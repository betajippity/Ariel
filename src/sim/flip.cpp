// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: flip.cpp
// Implements the FLIP sim

#include "flip.hpp"
#include "../math/kernels.inl"
#include <omp.h>

using namespace fluidCore;

flipsim::flipsim(const vec3& maxres, const float& density){
	dimensions = maxres;	
	pgrid = new particlegrid(maxres);
	max_density = 0.0f;
	this->density = density;


	float h = density/dimensions.x;
	FOR_EACH_CELL(10, 10, 10){
		particle* p = new particle;
		p->p = (vec3(i,j,k) + vec3(0.5f))*h;
		p->type = FLUID;
		p->mass = 1.0f;
		particles.push_back(p);
	}
	pgrid->sort(particles);
	// max_dens = 1.0;
}

flipsim::~flipsim(){
	delete pgrid;
	int particlecount = particles.size();
	for(int i=0; i<particlecount; i++){
		delete particles[i];
	}
	particles.clear();
	clearMacgrid(mgrid);
}

void flipsim::computeDensity(){

	int particlecount = particles.size();

	#pragma omp parallel for
	for(int i=0; i<particlecount; i++){
		//Find neighbours
		if(particles[i]->type==SOLID){
			particles[i]->density = 1.0f;
		}else{
			vec3 position = particles[i]->p;
			position = glm::max(vec3(0), glm::min(dimensions-vec3(1), dimensions*position));
			vector<particle *> neighbors = pgrid->getCellNeighbors(position, vec3(1));

			float weightsum = 0.0f;
			int neighborscount = neighbors.size();
			for(int m=0; m<neighborscount; m++){
				if(neighbors[m]->type!=SOLID){
					float sqd = mathCore::sqrlength(neighbors[m]->p, position);
					//TODO: figure out a better density smooth approx than density/dimensions.x
					float weight = neighbors[m]->mass * mathCore::smooth(sqd, 4.0f*density/dimensions.x);
					float weightsum = weightsum + weight;
				}
			}
			particles[i]->density = weightsum/max_density;
		}
	}
}