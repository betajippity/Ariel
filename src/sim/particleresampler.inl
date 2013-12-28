// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: particleresampler.inl
// Breakout file for some particle resampling operations stuff

#ifndef PARTICLERESAMPLER_INL
#define PARTICLERESAMPLER_INL

#include "../grid/macgrid.inl"
#include "../grid/particlegrid.hpp"
#include "../grid/levelset.hpp"
#include "../utilities/utilities.h"
#include "../grid/gridutils.inl"
#include <omp.h>

using namespace std;
using namespace glm;

namespace fluidCore {
//====================================
// Struct and Function Declarations
//====================================

//Forward declarations for externed inlineable methods
extern inline void resampleParticles(particlegrid* pgrid, vector<particle*>& particles, const float& dt, 
									 const float& re, const vec3& dimensions);
inline vec3 resample(particlegrid* pgrid, const vec3& p, const vec3& u, float re, const vec3& dimensions);


//====================================
// Function Implementations
//====================================

void resampleParticles(particlegrid* pgrid, vector<particle*>& particles, const float& dt, const float& re, 
					   const vec3& dimensions){
	int nx = (int)dimensions.x; int ny = (int)dimensions.y; int nz = (int)dimensions.z;
	float maxd = glm::max(glm::max(nx, ny), nz);
	pgrid->sort(particles);

	float springforce = 50.0f;

	//use springs to temporarily displace particles
	#pragma omp parallel for
	for(int n0=0; n0<particles.size(); n0++){	
		if(particles[n0]->type==FLUID){
			particle* p = particles[n0];
			vec3 spring(0.0f, 0.0f, 0.0f);
			float x = glm::max(0.0f,glm::min((float)maxd,maxd*p->p.x));
			float y = glm::max(0.0f,glm::min((float)maxd,maxd*p->p.y));
			float z = glm::max(0.0f,glm::min((float)maxd,maxd*p->p.z));
			vector<particle*> neighbors = pgrid->getCellNeighbors(vec3(x,y,z),vec3(1));
			for(int n1=0; n1<neighbors.size(); n1++){
				particle* np = neighbors[n1];
				if(p!=np){
					float dist = length(p->p-np->p);
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

	#pragma omp parallel for
	for(int n=0; n<particles.size(); n++){
		if(particles[n]->type == FLUID){
			particle* p = particles[n];
			p->t2.x = p->u.x;
			p->t2.y = p->u.y;
			p->t2.z = p->u.z;
			p->t2 = resample(pgrid, p->t, p->t2, re, dimensions);
		}
	}

	#pragma omp parallel for
	for(int n=0; n<particles.size(); n++){
		if(particles[n]->type == FLUID){
			particle *p = particles[n];
			p->p = p->t;
			p->u = p->t2;
		}
	}
}

vec3 resample(particlegrid* pgrid, const vec3& p, const vec3& u, float re, const vec3& dimensions){
	int nx = (int)dimensions.x; int ny = (int)dimensions.y; int nz = (int)dimensions.z;
	float maxd = glm::max(glm::max(nx, ny), nz);

	float wsum = 0.0f;
	vec3 ru = vec3(0);

	float x = glm::max(0.0f,glm::min((float)maxd-1,maxd*p.x));
	float y = glm::max(0.0f,glm::min((float)maxd-1,maxd*p.y));
	float z = glm::max(0.0f,glm::min((float)maxd-1,maxd*p.z));
	vector<particle*> neighbors = pgrid->getCellNeighbors(vec3(x,y,z),vec3(1));

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
