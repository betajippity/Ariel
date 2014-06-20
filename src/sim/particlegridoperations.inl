// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: particlegridsplat.inl
// Breakout file for some particle-grid operations stuff

#ifndef PARTICLEGRIDOPERATIONS_INL
#define PARTICLEGRIDOPERATIONS_INL

#include <tbb/tbb.h>
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
extern inline void SplatParticlesToMACGrid(ParticleGrid* sgrid, std::vector<Particle*>& particles,
										   MacGrid* mgrid);
extern inline void SplatMACGridToParticles(std::vector<Particle*>& particles, MacGrid* mgrid);
extern inline void EnforceBoundaryVelocity(MacGrid* mgrid);
extern inline glm::vec3 InterpolateVelocity(glm::vec3 p, MacGrid* mgrid);
inline float CheckWall(Grid<int>* A, const int& x, const int& y, const int& z);
inline float Interpolate(Grid<float>* q, glm::vec3 p, glm::vec3 n);
	
//====================================
// Function Implementations
//====================================

float CheckWall(Grid<int>* A, const int& x, const int& y, const int& z){
	if(A->GetCell(x,y,z)==SOLID){ //inside wall
		return 1.0f;
	}else{
		return -1.0f;
	}
}

void EnforceBoundaryVelocity(MacGrid* mgrid){
	unsigned int x = (unsigned int)mgrid->m_dimensions.x;
	unsigned int y = (unsigned int)mgrid->m_dimensions.y;
	unsigned int z = (unsigned int)mgrid->m_dimensions.z;

	//for every x face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x+1),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
			  	for(unsigned int j = 0; j < y; ++j){ 
			    	for(unsigned int k = 0; k < z; ++k){
			    		if(i==0 || i==x){
			    			mgrid->m_u_x->SetCell(i,j,k, 0.0f);
			    		}
						if( i<x && i>0 && CheckWall(mgrid->m_A, i, j, k)*
							CheckWall(mgrid->m_A, i-1, j, k) < 0 ) {
							mgrid->m_u_x->SetCell(i,j,k, 0.0f);
						}
			    	}
			   	}
			} 
		}
	);
	//for every y face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
			  	for(unsigned int j = 0; j < y+1; ++j){ 
			    	for(unsigned int k = 0; k < z; ++k){
			    		if(j==0 || j==y){
			    			mgrid->m_u_y->SetCell(i,j,k, 0.0f);
			    		}
						if( j<y && j>0 && CheckWall(mgrid->m_A, i, j, k)*
										  CheckWall(mgrid->m_A, i, j-1, k) < 0 ) {
							mgrid->m_u_y->SetCell(i,j,k, 0.0f);
						}
			    	}
			   	}
			} 
		}
	);
	//for every z face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
			  	for(unsigned int j = 0; j < y; ++j){ 
			    	for(unsigned int k = 0; k < z+1; ++k){
			    		if(k==0 || k==z){
			    			mgrid->m_u_z->SetCell(i,j,k, 0.0f);
			    		}
						if( k<z && k>0 && CheckWall(mgrid->m_A, i, j, k)*
										  CheckWall(mgrid->m_A, i, j, k-1) < 0 ) {
							mgrid->m_u_z->SetCell(i,j,k, 0.0f);
						}
			    	}
			   	}
			} 
		}
	);
}

float Interpolate(Grid<float>* q, glm::vec3 p, glm::vec3 n){
	float x = glm::max(0.0f,glm::min(n.x,p.x));
	float y = glm::max(0.0f,glm::min(n.y,p.y));
	float z = glm::max(0.0f,glm::min(n.z,p.z));
	int i = glm::min(x,n.x-2);
	int j = glm::min(y,n.y-2);
	int k = glm::min(z,n.z-2);
	float term1 = ((i+1-x)*q->GetCell(i,j,k)+(x-i)*q->GetCell(i+1,j,k))*(j+1-y);
	float term2 = ((i+1-x)*q->GetCell(i,j+1,k)+(x-i)*q->GetCell(i+1,j+1,k))*(y-j);
	float term3 = ((i+1-x)*q->GetCell(i,j,k+1)+(x-i)*q->GetCell(i+1,j,k+1))*(j+1-y);
	float term4 = ((i+1-x)*q->GetCell(i,j+1,k+1)+(x-i)*q->GetCell(i+1,j+1,k+1))*(y-j);
	return (k+1-z)*(term1 + term2) + (z-k)*(term3 + term4);
}

glm::vec3 InterpolateVelocity(glm::vec3 p, MacGrid* mgrid){
	int x = (int)mgrid->m_dimensions.x; int y = (int)mgrid->m_dimensions.y; 
	int z = (int)mgrid->m_dimensions.z;
	float maxd = glm::max(glm::max(x,y),z);
	x = maxd; y = maxd; z = maxd;
	glm::vec3 u;
	u.x = Interpolate(mgrid->m_u_x, glm::vec3(x*p.x, y*p.y-0.5f, z*p.z-0.5f), 
					  glm::vec3(x+1, y, z));
	u.y = Interpolate(mgrid->m_u_y, glm::vec3(x*p.x-0.5f, y*p.y, z*p.z-0.5f), 
					  glm::vec3(x, y+1, z));
	u.z = Interpolate(mgrid->m_u_z, glm::vec3(x*p.x-0.5f, y*p.y-0.5f, z*p.z), 
					  glm::vec3(x, y, z+1));
	return u;
}

void SplatMACGridToParticles(std::vector<Particle*>& particles, MacGrid* mgrid){
	unsigned int particleCount = particles.size();
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particleCount),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				particles[i]->m_u = InterpolateVelocity(particles[i]->m_p, mgrid);
			}
		}
	);
}

void SplatParticlesToMACGrid(ParticleGrid* sgrid, std::vector<Particle*>& particles, 
							 MacGrid* mgrid){
	
	float RE = 1.4f; //sharpen kernel weight

	int x = (int)mgrid->m_dimensions.x; int y = (int)mgrid->m_dimensions.y; 
	int z = (int)mgrid->m_dimensions.z;
	float maxd = glm::max(glm::max(x,y),z);

	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x+1),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				for(unsigned int j = 0; j < y+1; ++j){
					for(unsigned int k = 0; k < z+1; ++k){
						std::vector<Particle*> neighbors;
						//Splat X direction
						if(j<y && k<z){
							glm::vec3 px = glm::vec3(i, j+0.5f, k+0.5f);
							float sumw = 0.0f;
							float sumx = 0.0f;
							neighbors = sgrid->GetWallNeighbors(glm::vec3(i,j,k), 
																glm::vec3(1,2,2));
							for(unsigned int n=0; n<neighbors.size(); n++){
								Particle* p = neighbors[n];
								if(p->m_type == FLUID){
									glm::vec3 pos;
									pos.x = glm::max(0.0f,glm::min(maxd,maxd*p->m_p.x));
									pos.y = glm::max(0.0f,glm::min(maxd,maxd*p->m_p.y));
									pos.z = glm::max(0.0f,glm::min(maxd,maxd*p->m_p.z));
									float w = p->m_mass * mathCore::Sharpen(
														mathCore::Sqrlength(pos,px),RE);
									sumx += w*p->m_u.x;
									sumw += w;
								}
							}
							float uxsum = 0.0f;
							if(sumw>0){ 
								uxsum = sumx/sumw;
							}
							mgrid->m_u_x->SetCell(i,j,k,uxsum);
						}
						neighbors.clear();

						//Splat Y direction
						if(i<x && k<z){
							glm::vec3 py = glm::vec3(i+0.5f, j, k+0.5f);
							float sumw = 0.0f;
							float sumy = 0.0f;
							neighbors = sgrid->GetWallNeighbors(glm::vec3(i,j,k), 
																glm::vec3(2,1,2));
							for(unsigned int n=0; n<neighbors.size(); n++){
								Particle* p = neighbors[n];
								if(p->m_type == FLUID){
									glm::vec3 pos;
									pos.x = glm::max(0.0f,glm::min(maxd,maxd*p->m_p.x));
									pos.y = glm::max(0.0f,glm::min(maxd,maxd*p->m_p.y));
									pos.z = glm::max(0.0f,glm::min(maxd,maxd*p->m_p.z));
									float w = p->m_mass * mathCore::Sharpen(
														mathCore::Sqrlength(pos,py),RE);
									sumy += w*p->m_u.y;
									sumw += w;
								}
							}
							float uysum = 0.0f;
							if(sumw>0){
								uysum = sumy/sumw;
							}
							mgrid->m_u_y->SetCell(i,j,k,uysum);
						}
						neighbors.clear();

						//Splat Z direction
						if(i<x && j<y){
							glm::vec3 pz = glm::vec3(i+0.5f, j+0.5f, k);
							float sumw = 0.0f;
							float sumz = 0.0f;
							neighbors = sgrid->GetWallNeighbors(glm::vec3(i,j,k), 
																glm::vec3(2,2,1));
							for(unsigned int n=0; n<neighbors.size(); n++){
								Particle* p = neighbors[n];
								if(p->m_type == FLUID){
									glm::vec3 pos;
									pos.x = glm::max(0.0f,glm::min(maxd,maxd*p->m_p.x));
									pos.y = glm::max(0.0f,glm::min(maxd,maxd*p->m_p.y));
									pos.z = glm::max(0.0f,glm::min(maxd,maxd*p->m_p.z));
									float w = p->m_mass * mathCore::Sharpen(
														mathCore::Sqrlength(pos,pz),RE);
									sumz += w*p->m_u.z;
									sumw += w;
								}
							}
							float uzsum = 0.0f;
							if(sumw>0){
								uzsum = sumz/sumw;
							}
							mgrid->m_u_z->SetCell(i,j,k,uzsum);
						}
						neighbors.clear();
					}
				}
			}
		}
	);
}

}

#endif
