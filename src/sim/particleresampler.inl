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
extern inline void resampleParticles(ParticleGrid* pgrid, std::vector<particle*>& particles, 
									 const float& dt, const float& re, 
									 const glm::vec3& dimensions);
inline glm::vec3 resample(ParticleGrid* pgrid, const glm::vec3& p, const glm::vec3& u, float re, 
					 	  const glm::vec3& dimensions);


//====================================
// Function Implementations
//====================================

void resampleParticles(ParticleGrid* pgrid, std::vector<particle*>& particles, const float& dt, 
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
				if(particles[n0]->type==FLUID){
					particle* p = particles[n0];
					glm::vec3 spring(0.0f, 0.0f, 0.0f);
					float x = glm::max(0.0f,glm::min((float)maxd,maxd*p->p.x));
					float y = glm::max(0.0f,glm::min((float)maxd,maxd*p->p.y));
					float z = glm::max(0.0f,glm::min((float)maxd,maxd*p->p.z));
					std::vector<particle*> neighbors = pgrid->GetCellNeighbors(glm::vec3(x,y,z),
																			   glm::vec3(1));
					unsigned int neighborsCount = neighbors.size();
					for(unsigned int n1=0; n1<neighborsCount; ++n1){
						particle* np = neighbors[n1];
						if(p!=np){
							float dist = glm::length(p->p-np->p);
							float w = springforce * np->mass * mathCore::smooth(dist*dist,re);
							if(dist > 0.1f*re){
								spring.x += w * (p->p.x-np->p.x) / dist * re;
								spring.y += w * (p->p.y-np->p.y) / dist * re;
								spring.z += w * (p->p.z-np->p.z) / dist * re;
							}else{
								if(np->type == FLUID){
									spring.x += 0.01f*re/dt*(rand()%101)/100.0f;
									spring.y += 0.01f*re/dt*(rand()%101)/100.0f;
									spring.z += 0.01f*re/dt*(rand()%101)/100.0f;
								}else{
									spring.x += 0.05f*re/dt*np->n.x;
									spring.y += 0.05f*re/dt*np->n.y;
									spring.z += 0.05f*re/dt*np->n.z;
								}
							}
						}
					}
					p->t.x = p->p.x + dt*spring.x;
					p->t.y = p->p.y + dt*spring.y;
					p->t.z = p->p.z + dt*spring.z;
				}
			}
		}
	);

	particleCount = particles.size();
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particleCount),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int n=r.begin(); n!=r.end(); ++n){	
				if(particles[n]->type == FLUID){
					particle* p = particles[n];
					p->t2.x = p->u.x;
					p->t2.y = p->u.y;
					p->t2.z = p->u.z;
					p->t2 = resample(pgrid, p->t, p->t2, re, dimensions);
				}
			}
		}
	);

	particleCount = particles.size();
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particleCount),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int n=r.begin(); n!=r.end(); ++n){	
				if(particles[n]->type == FLUID){
					particle *p = particles[n];
					p->p = p->t;
					p->u = p->t2;
				}
			}
		}
	);
}

glm::vec3 resample(ParticleGrid* pgrid, const glm::vec3& p, const glm::vec3& u, float re, 
				   const glm::vec3& dimensions){
	int nx = (int)dimensions.x; int ny = (int)dimensions.y; int nz = (int)dimensions.z;
	float maxd = glm::max(glm::max(nx, ny), nz);

	float wsum = 0.0f;
	glm::vec3 ru = glm::vec3(0);

	float x = glm::max(0.0f,glm::min((float)maxd-1,maxd*p.x));
	float y = glm::max(0.0f,glm::min((float)maxd-1,maxd*p.y));
	float z = glm::max(0.0f,glm::min((float)maxd-1,maxd*p.z));
	std::vector<particle*> neighbors = pgrid->GetCellNeighbors(glm::vec3(x,y,z),glm::vec3(1));

	for(int n=0; n<neighbors.size(); n++){
		particle *np = neighbors[n];
		if(np->type == FLUID){
			float dist2 = mathCore::sqrlength(p,np->p);
			float w = np->mass * mathCore::sharpen(dist2,re);
			ru += w * np->u;
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
