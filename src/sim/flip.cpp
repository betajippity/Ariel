// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: flip.cpp
// Implements the FLIP sim

#include "flip.hpp"
#include "../math/kernels.inl"
#include "particlegridoperations.inl"
#include "solver.inl"
#include <omp.h>

using namespace fluidCore;

flipsim::flipsim(const vec3& maxres, sceneCore::scene* s, const float& density){
	dimensions = maxres;	
	pgrid = new particlegrid(maxres);
	mgrid = createMacgrid(maxres);
	max_density = 0.0f;
	this->density = density;
	scene = s;
	timestep = 0;
	stepsize = 0.005f;
	subcell = 1;
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

void flipsim::init(){
	//We need to figure out maximum particle pressure, so we generate a bunch of temporary particles
	//inside of a known area, sort them back onto the underlying grid, and calculate the density
	float maxd = glm::max(glm::max(dimensions.x, dimensions.z), dimensions.y);
	float h = density/maxd;
 	//generate temp particles
	for(int i = 0; i < 10; i++){  				//FOR_EACH_CELL
		for(int j = 0; j < 10; j++){ 
			for(int k = 0; k < 10; k++){ 
				particle* p = new particle;
				p->p = (vec3(i,j,k) + vec3(0.5f))*h;
				p->type = FLUID;
				p->mass = 1.0f;
				particles.push_back(p);
			}
		}
	}
	pgrid->sort(particles);
	max_density = 1.0f;
	computeDensity(); 
	max_density = 0.0f;
	//sum densities across particles
	for( int n=0; n<particles.size(); n++ ) {
		particle *p = particles[n];
		max_density = fmax(max_density,p->density);
		delete p;
	}
	particles.clear();
	cout << "Maxdensity: " << max_density << endl;

	//Generate particles
	scene->generateParticles(particles, dimensions, density, pgrid);

	//Sort particles to grid cells and remove fluid particles that are inside walls
	pgrid->sort(particles);
	pgrid->markCellTypes(particles, mgrid.A, density);

}

void flipsim::step(){
	timestep++;
	cout << "=========================" << endl;
	cout << "Step: " << timestep << endl;
	//Compute density
	cout << "Sorting/computing density..." << endl;
	pgrid->sort(particles);
	computeDensity();
	//Add forces
	cout << "Applying internal forces..." << endl;
	applyExternalForces(); 
	//figure out what cell each particle goes in
	cout << "Splatting particles to MAC Grid..." << endl;
	splatParticlesToMACGrid(pgrid, particles, mgrid, dimensions);
	//build liquid level set
	cout << "Building liquid level set..." << endl;
	scene->rebuildLiquidLevelSet(particles);

	scene->getLiquidLevelSet()->writeVDBGridToFile("test.vdb");


	//set solid-liquid interface velocities to zero
	cout << "Enforcing boundary velocities..." << endl;
	enforceBoundaryVelocity(mgrid, scene->getSolidLevelSet(), dimensions);
	//projection step
	// cout << "Running project step..." << endl;
	// project();
}

void flipsim::project(){
	int x = (int)dimensions.x;
	int y = (int)dimensions.y;
	int z = (int)dimensions.z;

	float maxd = glm::max(glm::max(dimensions.x, dimensions.z), dimensions.y);
	float h = 1.0f/maxd; //cell width

	cout << "Computing divergence..." << endl;
	//compute divergence per cell
	//for now run single threaded, multithreaded seems to cause VDB write issues here
	// #pragma omp parallel for
	for(int i = 0; i < x; i++){
		for(int j = 0; j < y; j++){
			for(int k = 0; k < z; k++){
				if(isCellFluid(i,j,k)==true){ //actually probably don't need a fluid check
					float divergence = (mgrid.u_x->getCell(i+1, j, k) - mgrid.u_x->getCell(i, j, k) + 
									    mgrid.u_y->getCell(i, j+1, k) - mgrid.u_y->getCell(i, j, k) +
									    mgrid.u_z->getCell(i, j, k+1) - mgrid.u_z->getCell(i, j, k)) / h;
					mgrid.D->setCell(i,j,k,divergence);
				}
			}
		}
	}

	cout << "Running solver..." << endl;
	solve(scene->getLiquidLevelSet(), scene->getSolidLevelSet(), mgrid, dimensions, subcell);

	//TODO: rest of project
}

void flipsim::applyExternalForces(){
	vec3 gravity = vec3(0,-9.8f, 0); //for now, just gravity
	int particlecount = particles.size();
	for(int i=0; i<particlecount; i++){
		particles[i]->u = gravity*stepsize;
	}
}

void flipsim::computeDensity(){

	float maxd = glm::max(glm::max(dimensions.x, dimensions.z), dimensions.y);

	int particlecount = particles.size();
	#pragma omp parallel for
	for(int i=0; i<particlecount; i++){
		//Find neighbours
		if(particles[i]->type==SOLID){
			particles[i]->density = 1.0f;
		}else{
			vec3 position = particles[i]->p;

			position.x = (int)fmax(0,fmin((int)dimensions.x-1,(int)dimensions.x*position.x));
			position.y = (int)fmax(0,fmin((int)dimensions.y-1,(int)dimensions.y*position.y));
			position.z = (int)fmax(0,fmin((int)dimensions.z-1,(int)dimensions.z*position.z));
			vector<particle *> neighbors;
			neighbors = pgrid->getCellNeighbors(position, vec3(1));
			float weightsum = 0.0f;
			int neighborscount = neighbors.size();
			for(int m=0; m<neighborscount; m++){
				if(neighbors[m]->type!=SOLID){
					float sqd = mathCore::sqrlength(neighbors[m]->p, particles[i]->p);
					//TODO: figure out a better density smooth approx than density/maxd
					float weight = neighbors[m]->mass * mathCore::smooth(sqd, 4.0f*density/maxd);
					weightsum = weightsum + weight;
				}
			}
			particles[i]->density = weightsum/max_density;
		}
	}
}

bool flipsim::isCellFluid(const int& x, const int& y, const int& z){
	if(scene->getLiquidLevelSet()->getCell(x,y,z)<0.0f &&
	   scene->getSolidLevelSet()->getCell(x,y,z)>=0.0f){
		return true;
	}else{
		return false;
	}
}

vector<particle*>* flipsim::getParticles(){
	return &particles;
}

vec3 flipsim::getDimensions(){
	return dimensions;
}

sceneCore::scene* flipsim::getScene(){
	return scene;	
}
