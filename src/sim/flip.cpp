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

	//Generate particles and sort
	scene->generateParticles(particles, dimensions, density, pgrid);
	pgrid->sort(particles);
	pgrid->markCellTypes(particles, mgrid.A, density);

	//Remove fluid particles that are stuck in walls
	for(vector<particle *>::iterator iter=particles.begin(); iter!=particles.end();) {
		particle &p = **iter;
		if(p.type == SOLID){
			iter++;
			continue;
		}
		int i = glm::min(dimensions.x-1,glm::max(0.0f,p.p.x*dimensions.x));
		int j = glm::min(dimensions.y-1,glm::max(0.0f,p.p.y*dimensions.y));
		int k = glm::min(dimensions.z-1,glm::max(0.0f,p.p.z*dimensions.z));
		if( mgrid.A->getCell(i,j,k) == SOLID ) {
			delete *iter;
			iter = particles.erase(iter);
		} else {
			iter ++;
		}
	}
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
	cout << "Applying external forces..." << endl;
	applyExternalForces(); 
	//figure out what cell each particle goes in
	cout << "Splatting particles to MAC Grid..." << endl;
	splatParticlesToMACGrid(pgrid, particles, mgrid);
	//set solid-liquid interface velocities to zero
	cout << "Enforcing boundary velocities..." << endl;
	enforceBoundaryVelocity(mgrid);
	//projection step
	cout << "Running project step..." << endl;
	project();
	cout << "Enforcing boundary velocities..." << endl;
	enforceBoundaryVelocity(mgrid);
	cout << "Extrapolating velocities..." << endl;
	extrapolateVelocity();
}

void flipsim::project(){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;

	float maxd = glm::max(glm::max(dimensions.x, dimensions.z), dimensions.y);
	float h = 1.0f/maxd; //cell width

	cout << "Computing divergence..." << endl;
	//compute divergence per cell
	// #pragma omp parallel for
	for(int i = 0; i < x; i++){
		for(int j = 0; j < y; j++){
			for(int k = 0; k < z; k++){
				float divergence = (mgrid.u_x->getCell(i+1, j, k) - mgrid.u_x->getCell(i, j, k) + 
								    mgrid.u_y->getCell(i, j+1, k) - mgrid.u_y->getCell(i, j, k) +
								    mgrid.u_z->getCell(i, j, k+1) - mgrid.u_z->getCell(i, j, k)) / h;
				mgrid.D->setCell(i,j,k,divergence);
			}
		}
	}

	cout << "Building liquid SDF..." << endl;
	//compute internal level set for liquid surface
	pgrid->buildSDF(mgrid, density);

	cout << "Running solver..." << endl;
	solve(mgrid, subcell);

	//subtract pressure gradient
	cout << "Subtracting pressure gradient..." << endl;
	subtractPressureGradient();
}

void flipsim::extrapolateVelocity(){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;

	floatgrid mark[3] = {floatgrid(), floatgrid(), floatgrid()};
	floatgrid wallmark[3] = {floatgrid(), floatgrid(), floatgrid()};

	//initalize temp grids with values
	//for every x face
	for(int i = 0; i < x+1; i++){ 	
	  	for(int j = 0; j < y; j++){ 
	    	for(int k = 0; k < z; k++){
				mark[0].setCell(i,j,k, (i>0 && mgrid.A->getCell(i-1,j,k)==FLUID) || 
									   (i<x && mgrid.A->getCell(i,j,k)==FLUID));
				wallmark[0].setCell(i,j,k, (i<=0 || mgrid.A->getCell(i-1,j,k)==SOLID) && 
										   (i>=x || mgrid.A->getCell(i,j,k)==SOLID));
	    	}
	    }
	}
	//for every y face
	for(int i = 0; i < x; i++){ 	
	  	for(int j = 0; j < y+1; j++){ 
	    	for(int k = 0; k < z; k++){
				mark[1].setCell(i,j,k, (j>0 && mgrid.A->getCell(i,j-1,k)==FLUID) || 
									   (j<y && mgrid.A->getCell(i,j,k)==FLUID));
				wallmark[1].setCell(i,j,k, (j<=0 || mgrid.A->getCell(i,j-1,k)==SOLID) && 
										   (j>=y || mgrid.A->getCell(i,j,k)==SOLID));
	    	}
	    }
	}
	//for every z face
	for(int i = 0; i < x; i++){ 	
	  	for(int j = 0; j < y; j++){ 
	    	for(int k = 0; k < z+1; k++){
				mark[2].setCell(i,j,k, (k>0 && mgrid.A->getCell(i,j,k-1)==FLUID) || 
									   (k<z && mgrid.A->getCell(i,j,k)==FLUID));
				wallmark[2].setCell(i,j,k, (k<=0 || mgrid.A->getCell(i,j,k-1)==SOLID) && 
										   (k>=z || mgrid.A->getCell(i,j,k)==SOLID));
	    	}
	    }
	}

	//extrapolate
	for(int i = 0; i < x+1; i++){ 	
	  	for(int j = 0; j < y+1; j++){ 
	    	for(int k = 0; k < z+1; k++){
				for(int n=0; n<3; n++){
					if(n!=0 && i>x-1){ continue; };
					if(n!=1 && j>y-1){ continue; };
					if(n!=2 && k>z-1){ continue; };
					if(!mark[n].getCell(i,j,k) && wallmark[n].getCell(i,j,k)){
						int wsum = 0;
						float sum = 0.0f;
						vec3 q[6] = { vec3(i-1,j,k), vec3(i+1,j,k), vec3(i,j-1,k), 
									  vec3(i,j+1,k), vec3(i,j,k-1), vec3(i,j,k+1) };
						for(int qk=0; qk<6; qk++){
							if(q[qk][0]>=0 && q[qk][0]<x+(n==0) && q[qk][1]>=0 && 
							   q[qk][1]<y+(n==1) && q[qk][2]>=0 && q[qk][2]<z+(n==2) ) {
								if(mark[n].getCell(q[qk][0],q[qk][1],q[qk][2])){
									wsum ++;
									if(n==0){
										sum += mgrid.u_x->getCell(q[qk][0],q[qk][1],q[qk][2]);
									}else if(n==1){
										sum += mgrid.u_y->getCell(q[qk][0],q[qk][1],q[qk][2]);
									}else if(n==2){
										sum += mgrid.u_z->getCell(q[qk][0],q[qk][1],q[qk][2]);
									}
								}
							}
						}
						if(wsum){
							if(n==0){
								mgrid.u_x->setCell(i,j,k,sum/wsum);
							}else if(n==1){
								mgrid.u_y->setCell(i,j,k,sum/wsum);
							}else if(n==2){
								mgrid.u_z->setCell(i,j,k,sum/wsum);
							}
						}
					}
				}
	    	}
	    }
	}
}

void flipsim::subtractPressureGradient(){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;

	float maxd = glm::max(glm::max(dimensions.x, dimensions.z), dimensions.y);
	float h = 1.0f/maxd; //cell width

	//for every x face
	#pragma omp parallel for
	for(int i = 0; i < x+1; i++){  
	  	for(int j = 0; j < y; j++){ 
	    	for(int k = 0; k < z; k++){
	    		if(i>0 && i<x){
					float pf = mgrid.P->getCell(i,j,k);
					float pb = mgrid.P->getCell(i-1,j,k);
					if(subcell && mgrid.L->getCell(i,j,k) * mgrid.L->getCell(i-1,j,k) < 0.0f){
						if(mgrid.L->getCell(i,j,k)<0.0f){
							pf = mgrid.P->getCell(i,j,k);
						}else{
							pf = mgrid.L->getCell(i,j,k)/
								 glm::min(1.0e-3f,mgrid.L->getCell(i-1,j,k))*mgrid.P->getCell(i-1,j,k);
						}
						if(mgrid.L->getCell(i-1,j,k)<0.0f){
							pb = mgrid.P->getCell(i-1,j,k);
						}else{
							pb = mgrid.L->getCell(i-1,j,k)/
								 glm::min(1.0e-6f,mgrid.L->getCell(i,j,k))*mgrid.P->getCell(i,j,k);
						}				
					}
					float xval = mgrid.u_x->getCell(i,j,k);
					xval -= (pf-pb)/h;
					mgrid.u_x->setCell(i,j,k,xval);
				}
	    	}
	    }
	} 
	//for every y face
	#pragma omp parallel for
 	for(int i = 0; i < x; i++){
	  	for(int j = 0; j < y+1; j++){
	    	for(int k = 0; k < z; k++){
	    		if(j>0 && j<y){
					float pf = mgrid.P->getCell(i,j,k);
					float pb = mgrid.P->getCell(i,j-1,k);   
					if(subcell && mgrid.L->getCell(i,j,k) * mgrid.L->getCell(i,j-1,k) < 0.0f){
						if(mgrid.L->getCell(i,j,k)<0.0f){
							pf = mgrid.P->getCell(i,j,k);
						}else{
							pf = mgrid.L->getCell(i,j,k)/
								 glm::min(1.0e-3f,mgrid.L->getCell(i,j-1,k))*mgrid.P->getCell(i,j-1,k);
						}
						if(mgrid.L->getCell(i,j-1,k)<0.0f){
							pb = mgrid.P->getCell(i,j-1,k);
						}else{
							pb = mgrid.L->getCell(i,j-1,k)/
								 glm::min(1.0e-6f,mgrid.L->getCell(i,j,k))*mgrid.P->getCell(i,j,k);
						}	
					}
					float yval = mgrid.u_y->getCell(i,j,k);
					yval -= (pf-pb)/h;
					mgrid.u_y->setCell(i,j,k,yval);
	    		} 		
	    	} 
	    }
	}
	//for every z face
	#pragma omp parallel for
 	for(int i = 0; i < x; i++){
	  	for(int j = 0; j < y; j++){
	    	for(int k = 0; k < z+1; k++){
	    		if(k>0 && k<z){
					float pf = mgrid.P->getCell(i,j,k);
					float pb = mgrid.P->getCell(i,j,k-1);
					if(subcell && mgrid.L->getCell(i,j,k) * mgrid.L->getCell(i,j,k-1) < 0.0f){
						if(mgrid.L->getCell(i,j,k)<0.0f){
							pf = mgrid.P->getCell(i,j,k);
						}else{
							pf = mgrid.L->getCell(i,j,k)/
								 glm::min(1.0e-3f,mgrid.L->getCell(i,j,k-1))*mgrid.P->getCell(i,j,k-1);
						}
						if(mgrid.L->getCell(i,j,k-1)<0.0f){
							pb = mgrid.P->getCell(i,j,k-1);
						}else{
							pb = mgrid.L->getCell(i,j,k-1)/
								 glm::min(1.0e-6f,mgrid.L->getCell(i,j,k))*mgrid.P->getCell(i,j,k);
						}	
					}
					float zval = mgrid.u_z->getCell(i,j,k);
					zval -= (pf-pb)/h;
					mgrid.u_z->setCell(i,j,k,zval);
	    		} 		
	    	} 
	    }
	}
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
