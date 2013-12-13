// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: particlegridsplat.inl
// Breakout file for some particle-grid splatting/mapping stuff

#ifndef PARTICLEGRIDSPLAT_INL
#define PARTICLEGRIDSPLAT_INL

#include "../grid/macgrid.inl"
#include "../grid/particlegrid.hpp"
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
										   macgrid& mgrid, const vec3& dimensions);
	
//====================================
// Function Implementations
//====================================

void splatParticlesToMACGrid(particlegrid* sgrid, vector<particle*>& particles, macgrid& mgrid, 
							 const vec3& dimensions){
	
	float RE = 1.4f; //sharpen kernel weight

	// float xtotal, ytotal, ztotal;
	// xtotal = 0;
	// ytotal = 0;
	// ztotal = 0;

	#pragma omp parallel for
	for(int i = 0; i < (int)dimensions.x+1; i++){
		for(int j = 0; j < (int)dimensions.y+1; j++){
			for(int k = 0; k < (int)dimensions.z+1; k++){
				vector<particle*> neighbors;
				//Splat X direction
				if(j<dimensions.y && k<dimensions.z){
					vec3 px = vec3(i, j+0.5f, k+0.5f);
					float sumw = 0.0f;
					float sumx = 0.0f;
					neighbors = sgrid->getWallNeighbors(vec3(i,j,k), vec3(1,2,2));
					for(int n=0; n<neighbors.size(); n++){
						particle* p = neighbors[n];
						vec3 pos;
						pos.x = fmax(0,fmin((int)dimensions.x,(int)dimensions.x*p->p.x));
						pos.y = fmax(0,fmin((int)dimensions.y,(int)dimensions.y*p->p.y));
						pos.z = fmax(0,fmin((int)dimensions.z,(int)dimensions.z*p->p.z));
						float w = p->mass * mathCore::sharpen(mathCore::sqrlength(pos,px),RE);
						sumx += w*p->u.x;
						sumw += w;
					}
					float uxsum = 0.0f;
					if(sumw>0){ //hope this works
						uxsum = sumx/sumw;
					}
					mgrid.u_x->setCell(i,j,k,uxsum);
					// xtotal = xtotal + uxsum;
				}
				neighbors.clear();

				//Splat Y direction
				if(i<dimensions.x && k<dimensions.z){
					vec3 py = vec3(i+0.5f, j, k+0.5f);
					float sumw = 0.0f;
					float sumy = 0.0f;
					neighbors = sgrid->getWallNeighbors(vec3(i,j,k), vec3(2,1,2));
					for(int n=0; n<neighbors.size(); n++){
						particle* p = neighbors[n];
						vec3 pos;
						pos.x = fmax(0,fmin((int)dimensions.x,(int)dimensions.x*p->p.x));
						pos.y = fmax(0,fmin((int)dimensions.y,(int)dimensions.y*p->p.y));
						pos.z = fmax(0,fmin((int)dimensions.z,(int)dimensions.z*p->p.z));
						float w = p->mass * mathCore::sharpen(mathCore::sqrlength(pos,py),RE);
						sumy += w*p->u.y;
						sumw += w;
					}
					float uysum = 0.0f;
					if(sumw>0){ //hope this works
						uysum = sumy/sumw;
					}
					mgrid.u_y->setCell(i,j,k,uysum);
					// ytotal = ytotal + uysum;
				}
				neighbors.clear();

				//Splat Z direction
				if(i<dimensions.x && j<dimensions.y){
					vec3 pz = vec3(i+0.5f, j+0.5f, k);
					float sumw = 0.0f;
					float sumz = 0.0f;
					neighbors = sgrid->getWallNeighbors(vec3(i,j,k), vec3(2,2,1));
					for(int n=0; n<neighbors.size(); n++){
						particle* p = neighbors[n];
						vec3 pos;
						pos.x = fmax(0,fmin((int)dimensions.x,(int)dimensions.x*p->p.x));
						pos.y = fmax(0,fmin((int)dimensions.y,(int)dimensions.y*p->p.y));
						pos.z = fmax(0,fmin((int)dimensions.z,(int)dimensions.z*p->p.z));
						float w = p->mass * mathCore::sharpen(mathCore::sqrlength(pos,pz),RE);
						sumz += w*p->u.z;
						sumw += w;
					}
					float uzsum = 0.0f;
					if(sumw>0){ //hope this works
						uzsum = sumz/sumw;
					}
					mgrid.u_z->setCell(i,j,k,uzsum);
					// ztotal = ztotal + uzsum;
				}
				neighbors.clear();
			}
		}
	}
	// cout << xtotal << " " << ytotal << " " << ztotal << endl;
}

}

#endif
