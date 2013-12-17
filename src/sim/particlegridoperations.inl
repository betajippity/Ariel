// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: particlegridsplat.inl
// Breakout file for some particle-grid operations stuff

#ifndef PARTICLEGRIDOPERATIONS_INL
#define PARTICLEGRIDOPERATIONS_INL

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
extern inline void splatParticlesToMACGrid(particlegrid* sgrid, vector<particle*>& particles,
										   macgrid& mgrid);
extern inline void enforceBoundaryVelocity(macgrid& mgrid);
extern inline float checkWall(intgrid* A, const int& x, const int& y, const int& z);
	
//====================================
// Function Implementations
//====================================

float checkWall(intgrid* A, const int& x, const int& y, const int& z){
	if(A->getCell(x,y,z)==SOLID){ //inside wall
		return 1.0f;
	}else{
		return -1.0f;
	}
}

void enforceBoundaryVelocity(macgrid& mgrid){
	int x = (int)mgrid.dimensions.x;
	int y = (int)mgrid.dimensions.y;
	int z = (int)mgrid.dimensions.z;

	//for every x face
	#pragma omp parallel for
	for(int i = 0; i < x+1; i++){  
	  	for(int j = 0; j < y; j++){ 
	    	for(int k = 0; k < z; k++){
	    		if(i==0 || i==x){
	    			mgrid.u_x->setCell(i,j,k, 0.0f);
	    		}
				if( i<x && i>0 && checkWall(mgrid.A, i, j, k)*checkWall(mgrid.A, i-1, j, k) < 0 ) {
					mgrid.u_x->setCell(i,j,k, 0.0f);
				}
	    	}
	   	}
	} 
	//for every y face
	#pragma omp parallel for
	for(int i = 0; i < x; i++){  
	  	for(int j = 0; j < y+1; j++){ 
	    	for(int k = 0; k < z; k++){
	    		if(j==0 || j==y){
	    			mgrid.u_y->setCell(i,j,k, 0.0f);
	    		}
				if( j<y && j>0 && checkWall(mgrid.A, i, j, k)*checkWall(mgrid.A, i, j-1, k) < 0 ) {
					mgrid.u_y->setCell(i,j,k, 0.0f);
				}
	    	}
	   	}
	} 
	//for every z face
	#pragma omp parallel for
	for(int i = 0; i < x; i++){  
	  	for(int j = 0; j < y; j++){ 
	    	for(int k = 0; k < z+1; k++){
	    		if(k==0 || k==z){
	    			mgrid.u_z->setCell(i,j,k, 0.0f);
	    		}
				if( k<z && k>0 && checkWall(mgrid.A, i, j, k)*checkWall(mgrid.A, i, j, k-1) < 0 ) {
					mgrid.u_z->setCell(i,j,k, 0.0f);
				}
	    	}
	   	}
	} 
}

void splatParticlesToMACGrid(particlegrid* sgrid, vector<particle*>& particles, macgrid& mgrid){
	
	float RE = 1.4f; //sharpen kernel weight

	int x = (int)mgrid.dimensions.x; int y = (int)mgrid.dimensions.y; int z = (int)mgrid.dimensions.z;
	
	// openmp disabled here since VDB likes to throw an assert fail here when multithreaded for some reason
	// #pragma omp parallel for
	for(int i = 0; i < x+1; i++){
		for(int j = 0; j < y+1; j++){
			for(int k = 0; k < z+1; k++){
				vector<particle*> neighbors;
				//Splat X direction
				if(j<y && k<z){
					vec3 px = vec3(i, j+0.5f, k+0.5f);
					float sumw = 0.0f;
					float sumx = 0.0f;
					neighbors = sgrid->getWallNeighbors(vec3(i,j,k), vec3(1,2,2));
					for(int n=0; n<neighbors.size(); n++){
						particle* p = neighbors[n];
						if(p->type == FLUID){
							vec3 pos;
							pos.x = fmax(0,fmin(x,x*p->p.x));
							pos.y = fmax(0,fmin(y,y*p->p.y));
							pos.z = fmax(0,fmin(z,z*p->p.z));
							float w = p->mass * mathCore::sharpen(mathCore::sqrlength(pos,px),RE);
							sumx += w*p->u.x;
							sumw += w;
						}
					}
					float uxsum = 0.0f;
					if(sumw>0){ 
						uxsum = sumx/sumw;
					}
					mgrid.u_x->setCell(i,j,k,uxsum);
				}
				neighbors.clear();

				//Splat Y direction
				if(i<x && k<z){
					vec3 py = vec3(i+0.5f, j, k+0.5f);
					float sumw = 0.0f;
					float sumy = 0.0f;
					neighbors = sgrid->getWallNeighbors(vec3(i,j,k), vec3(2,1,2));
					for(int n=0; n<neighbors.size(); n++){
						particle* p = neighbors[n];
						if(p->type == FLUID){
							vec3 pos;
							pos.x = fmax(0,fmin(x,x*p->p.x));
							pos.y = fmax(0,fmin(y,y*p->p.y));
							pos.z = fmax(0,fmin(z,z*p->p.z));
							float w = p->mass * mathCore::sharpen(mathCore::sqrlength(pos,py),RE);
							sumy += w*p->u.y;
							sumw += w;
						}
					}
					float uysum = 0.0f;
					if(sumw>0){
						uysum = sumy/sumw;
					}
					mgrid.u_y->setCell(i,j,k,uysum);
				}
				neighbors.clear();

				//Splat Z direction
				if(i<x && j<y){
					vec3 pz = vec3(i+0.5f, j+0.5f, k);
					float sumw = 0.0f;
					float sumz = 0.0f;
					neighbors = sgrid->getWallNeighbors(vec3(i,j,k), vec3(2,2,1));
					for(int n=0; n<neighbors.size(); n++){
						particle* p = neighbors[n];
						if(p->type == FLUID){
							vec3 pos;
							pos.x = fmax(0,fmin(x,x*p->p.x));
							pos.y = fmax(0,fmin(y,y*p->p.y));
							pos.z = fmax(0,fmin(z,z*p->p.z));
							float w = p->mass * mathCore::sharpen(mathCore::sqrlength(pos,pz),RE);
							sumz += w*p->u.z;
							sumw += w;
						}
					}
					float uzsum = 0.0f;
					if(sumw>0){
						uzsum = sumz/sumw;
					}
					mgrid.u_z->setCell(i,j,k,uzsum);
				}
				neighbors.clear();
			}
		}
	}
}

}

#endif
