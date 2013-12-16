// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: solver.inl
// Breakout file for PCG solver

#ifndef SOLVER_INL
#define SOLVER_INL

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
extern inline void solve(levelset* fls, levelset* sls, macgrid& mgrid, const vec3& dimensions, 
						 const int& subcell);
extern inline void flipDivergence(macgrid& mgrid, const vec3& dimensions);
extern inline void buildPreconditioner(floatgrid* pc, levelset* L, levelset* S, const vec3& dimensions, 
									   const int& subcell);
extern inline float fluidRef(levelset* L, int i, int j, int k, int qi, int qj, int qk, vec3 dimensions);
extern inline float preconditionerRef(floatgrid* p, int i, int j, int k, vec3 dimensions);
extern inline float fluidDiag(levelset* L, levelset* S, int i, int j, int k, vec3 dimensions, int subcell);
//====================================
// Function Implementations
//====================================

//Helper for preconditioner builder
float fluidRef(levelset* L, int i, int j, int k, int qi, int qj, int qk, vec3 dimensions){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	if( i<0 || i>x-1 || j<0 || j>y-1 || k<0 || k>z-1 || L->getCell(i,j,k)>=0.0f ){ //if not liquid
		return 0.0;
	} 
	if( qi<0 || qi>x-1 || qj<0 || qj>y-1 || qk<0 || qk>z-1 || L->getCell(i,j,k)>=0.0f ){ //if not liquid
		return 0.0;
	} 
	return -1.0;	
}

//Helper for preconditioner builder
float preconditionerRef(floatgrid* p, int i, int j, int k, vec3 dimensions){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	if( i<0 || i>x-1 || j<0 || j>y-1 || k<0 || k>z-1 || p->getCell(i,j,k)>=0.0f ){ //if not liquid
		return 0.0f;
	} 
	return p->getCell(i,j,k);
}

//Helper for preconditioner builder
float fluidDiag(levelset* L, levelset* S, int i, int j, int k, vec3 dimensions, int subcell){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	float diag = 6.0f;
	if( L->getCell(i,j,k) >= 0.0f ){ //if not liquid
		return diag;
	}
	vec3 q[] = { vec3(i-1,j,k), vec3(i+1,j,k), vec3(i,j-1,k), 
				 vec3(i,j+1,k), vec3(i,j,k-1), vec3(i,j,k+1) };
	for(int m=0; m<6; m++){
		int qi = q[m][0];
		int qj = q[m][1];
		int qk = q[m][2];
		if( qi<0 || qi>x-1 || qj<0 || qj>y-1 || qk<0 || qk>z-1 || S->getCell(i,j,k)<0.0f ){ //if in wall
			diag -= 1.0f;
			cout << "a " << diag << endl;
		}else if(L->getCell(i,j,k) >= 0.0f && subcell){ //if not liquid aka if air
			diag -= L->getCell(qi,qj,qk)/glm::min(1.0e-6f,L->getCell(i,j,k));
			cout << "b " << diag << endl;
		}
	}
	cout << diag << endl;
	return diag;
}

void buildPreconditioner(floatgrid* pc, levelset* L, levelset* S, const vec3& dimensions, 
						 const int& subcell){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	float a = 0.25f;
	//for now run single threaded, multithreaded seems to cause VDB write issues here
	// #pragma omp parallel for
	for( int gn=0; gn<x*y*z; gn++ ) { 
		int i=(gn%((x)*(y)))%(z); int j=(gn%((x)*(y)))/(z); int k = gn/((x)*(y)); 

		if(L->getCell(i,j,k)<0.0f){			//If fluid
			float left = fluidRef(L,i-1,j,k,i,j,k,dimensions)*preconditionerRef(pc,i-1,j,k,dimensions);
			float bottom = fluidRef(L,i,j-1,k,i,j,k,dimensions)*preconditionerRef(pc,i,j-1,k,dimensions);
			float back = fluidRef(L,i,j,k-1,i,j,k,dimensions)*preconditionerRef(pc,i,j,k-1,dimensions);
			float diag = fluidDiag(L,S,i,j,k,dimensions, subcell);
			double e = diag - (left*left) - (bottom*bottom) - (back*back);
			if( e < a*diag ){
				e = diag;
			}
			pc->setCell(i,j,k, 1.0f/sqrt(e));
		}
	}
}

void flipDivergence(macgrid& mgrid, const vec3& dimensions){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	#pragma omp parallel for
	for( int gn=0; gn<x*y*z; gn++ ) { 
		int i=(gn%((x)*(y)))%(z); int j=(gn%((x)*(y)))/(z); int k = gn/((x)*(y)); 
		float flippedD = -mgrid.D->getCell(i,j,k);
		mgrid.D->setCell(i,j,k,flippedD);
	}
}

void solve(levelset* fls, levelset* sls, macgrid& mgrid, const vec3& dimensions, const int& subcell){
	//flip divergence
	cout << "Flipping divergence..." << endl;
	flipDivergence(mgrid, dimensions);

	//build preconditioner
	cout << "Building preconditioner matrix..." << endl;
	floatgrid* preconditioner = new floatgrid(0.0f);
	buildPreconditioner(preconditioner, fls, sls, dimensions, subcell);

	delete preconditioner;
}

}

#endif
