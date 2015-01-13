// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: particleresampler.inl
// Breakout file for some particle resampling operations stuff

#ifndef PARTICLERESAMPLER_INL
#define PARTICLERESAMPLER_INL

#include <tbb/tbb.h>,
#include "../grid/macgrid.inl"
#include "../grid/particlegrid.hpp"
#include "../grid/levelset.hpp"
#include "../utilities/utilities.h"
#include "../grid/gridutils.inl"

namespace fluidCore {
//====================================
// Struct and Function Declarations
//====================================

//Forward declarations for externed inlineable methods
extern inline void ResampleParticles(ParticleGrid* pgrid, 
                                     std::vector<Particle*>& particles, 
                                     sceneCore::Scene* scene, const float& frame, const float& dt,
                                     const float& re, const glm::vec3& dimensions);
inline glm::vec3 Resample(ParticleGrid* pgrid, const glm::vec3& p, const glm::vec3& u, float re, 
                          const glm::vec3& dimensions);


//====================================
// Function Implementations
//====================================

void ResampleParticles(ParticleGrid* pgrid, std::vector<Particle*>& particles,
                       sceneCore::Scene* scene, const float& frame, const float& dt, 
                       const float& re, const glm::vec3& dimensions){
    int nx = (int)dimensions.x; int ny = (int)dimensions.y; int nz = (int)dimensions.z;
    float maxd = glm::max(glm::max(nx, ny), nz);
    pgrid->Sort(particles);

    float springforce = 50.0f;

    //use springs to temporarily displace particles
    unsigned int particleCount = particles.size();
    tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particleCount),
        [=](const tbb::blocked_range<unsigned int>& r){
            for(unsigned int n0=r.begin(); n0!=r.end(); ++n0){  
                if(particles[n0]->m_type==FLUID){
                    Particle* p = particles[n0];
                    glm::vec3 spring(0.0f, 0.0f, 0.0f);
                    float x = glm::max(0.0f,glm::min((float)maxd,maxd*p->m_p.x));
                    float y = glm::max(0.0f,glm::min((float)maxd,maxd*p->m_p.y));
                    float z = glm::max(0.0f,glm::min((float)maxd,maxd*p->m_p.z));
                    std::vector<Particle*> neighbors = pgrid->GetCellNeighbors(glm::vec3(x,y,z),
                                                                               glm::vec3(1));
                    unsigned int neighborsCount = neighbors.size();
                    for(unsigned int n1=0; n1<neighborsCount; ++n1){
                        Particle* np = neighbors[n1];
                        if(p!=np){
                            float dist = glm::length(p->m_p-np->m_p);
                            float w = springforce * np->m_mass * mathCore::Smooth(dist*dist,re);
                            if(dist > 0.1f*re){
                                spring.x += w * (p->m_p.x-np->m_p.x) / dist * re;
                                spring.y += w * (p->m_p.y-np->m_p.y) / dist * re;
                                spring.z += w * (p->m_p.z-np->m_p.z) / dist * re;
                            }else{
                                if(np->m_type == FLUID){
                                    spring.x += 0.01f*re/dt*(rand()%101)/100.0f;
                                    spring.y += 0.01f*re/dt*(rand()%101)/100.0f;
                                    spring.z += 0.01f*re/dt*(rand()%101)/100.0f;
                                }else{
                                    spring.x += 0.05f*re/dt*np->m_n.x;
                                    spring.y += 0.05f*re/dt*np->m_n.y;
                                    spring.z += 0.05f*re/dt*np->m_n.z;
                                }
                            }
                        }
                    }
                    p->m_t.x = p->m_p.x + dt*spring.x;
                    p->m_t.y = p->m_p.y + dt*spring.y;
                    p->m_t.z = p->m_p.z + dt*spring.z;
                }
            }
        }
    );

    particleCount = particles.size();
    tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particleCount),
        [=](const tbb::blocked_range<unsigned int>& r){
            for(unsigned int n=r.begin(); n!=r.end(); ++n){ 
                if(particles[n]->m_type == FLUID){
                    Particle* p = particles[n];
                    p->m_t2.x = p->m_u.x;
                    p->m_t2.y = p->m_u.y;
                    p->m_t2.z = p->m_u.z;
                    p->m_t2 = Resample(pgrid, p->m_t, p->m_t2, re, dimensions);
                }
            }
        }
    );

    particleCount = particles.size();
    tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particleCount),
        [=](const tbb::blocked_range<unsigned int>& r){
            for(unsigned int n=r.begin(); n!=r.end(); ++n){ 
                if(particles[n]->m_type == FLUID){
                    Particle *p = particles[n];
                    unsigned int solidGeomID = 0;
                    if(scene->CheckPointInsideSolidGeom(p->m_t*maxd, frame, solidGeomID)==false){
                        p->m_p = p->m_t;
                        p->m_u = p->m_t2;
                    }
                }
            }
        }
    );
}

glm::vec3 Resample(ParticleGrid* pgrid, const glm::vec3& p, const glm::vec3& u, float re, 
                   const glm::vec3& dimensions){
    int nx = (int)dimensions.x; int ny = (int)dimensions.y; int nz = (int)dimensions.z;
    float maxd = glm::max(glm::max(nx, ny), nz);

    float wsum = 0.0f;
    glm::vec3 ru = glm::vec3(0);

    float x = glm::max(0.0f,glm::min((float)maxd-1,maxd*p.x));
    float y = glm::max(0.0f,glm::min((float)maxd-1,maxd*p.y));
    float z = glm::max(0.0f,glm::min((float)maxd-1,maxd*p.z));
    std::vector<Particle*> neighbors = pgrid->GetCellNeighbors(glm::vec3(x,y,z),glm::vec3(1));

    for(int n=0; n<neighbors.size(); n++){
        Particle *np = neighbors[n];
        if(np->m_type == FLUID){
            float dist2 = mathCore::Sqrlength(p,np->m_p);
            float w = np->m_mass * mathCore::Sharpen(dist2,re);
            ru += w * np->m_u;
            wsum += w;
        }
    }
    if(wsum){
        ru /= wsum;
    } else {
        ru = u;
    }
    return ru;
}
}

#endif
