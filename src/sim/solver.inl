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
extern inline void solve(macgrid& mgrid, const int& subcell);
extern inline void flipDivergence(macgrid& mgrid);
extern inline void buildPreconditioner(floatgrid* pc, macgrid& mgrid, const int& subcell);
extern inline float fluidRef(intgrid* A, int i, int j, int k, int qi, int qj, int qk, vec3 dimensions);
extern inline float preconditionerRef(floatgrid* p, int i, int j, int k, vec3 dimensions);
extern inline float fluidDiag(intgrid* A, floatgrid* L, int i, int j, int k, vec3 dimensions, int subcell);
//====================================
// Function Implementations
//====================================

//Helper for preconditioner builder
float fluidRef(intgrid* A, int i, int j, int k, int qi, int qj, int qk, vec3 dimensions){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	if( i<0 || i>x-1 || j<0 || j>y-1 || k<0 || k>z-1 || A->getCell(i,j,k)!=FLUID ){ //if not liquid
		return 0.0;
	} 
	if( qi<0 || qi>x-1 || qj<0 || qj>y-1 || qk<0 || qk>z-1 || A->getCell(qi,qj,qk)!=FLUID ){ //if not liquid
		return 0.0;
	} 
	return -1.0;	
}

//Helper for preconditioner builder
float preconditionerRef(floatgrid* p, int i, int j, int k, vec3 dimensions){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	if( i<0 || i>x-1 || j<0 || j>y-1 || k<0 || k>z-1 || p->getCell(i,j,k)!=FLUID ){ //if not liquid
		return 0.0f;
	} 
	return p->getCell(i,j,k);
}

//Helper for preconditioner builder
float fluidDiag(intgrid* A, floatgrid* L, int i, int j, int k, vec3 dimensions, int subcell){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	float diag = 6.0;
	if( A->getCell(i,j,k) != FLUID ) return diag;
	int q[][3] = { {i-1,j,k}, {i+1,j,k}, {i,j-1,k}, {i,j+1,k}, {i,j,k-1}, {i,j,k+1} };
	for( int m=0; m<6; m++ ) {
		int qi = q[m][0];
		int qj = q[m][1];
		int qk = q[m][2];
		if( qi<0 || qi>x-1 || qj<0 || qj>y-1 || qk<0 || qk>z-1 || A->getCell(qi,qj,qk)==SOLID ){
			diag -= 1.0;
		}
		else if( A->getCell(qi,qj,qk)==AIR && subcell ) {
			diag -= L->getCell(qi,qj,qk)/glm::min(1.0e-6f,L->getCell(i,j,k));
		}
	}
	
	return diag;
}

void buildPreconditioner(floatgrid* pc, macgrid& mgrid, const int& subcell){
	int x = (int)mgrid.dimensions.x; int y = (int)mgrid.dimensions.y; int z = (int)mgrid.dimensions.z;
	float a = 0.25f;
	//for now run single threaded, multithreaded seems to cause VDB write issues here
	// #pragma omp parallel for
	for( int gn=0; gn<x*y*z; gn++ ) { 
		int i=(gn%((x)*(y)))%(z); int j=(gn%((x)*(y)))/(z); int k = gn/((x)*(y)); 
		if(mgrid.A->getCell(i,j,k)==FLUID){			//If fluid

			float left = fluidRef(mgrid.A,i-1,j,k,i,j,k,mgrid.dimensions)*
						 preconditionerRef(pc,i-1,j,k,mgrid.dimensions);
			float bottom = fluidRef(mgrid.A,i,j-1,k,i,j,k,mgrid.dimensions)*
						   preconditionerRef(pc,i,j-1,k,mgrid.dimensions);
			float back = fluidRef(mgrid.A,i,j,k-1,i,j,k,mgrid.dimensions)*
						 preconditionerRef(pc,i,j,k-1,mgrid.dimensions);
			float diag = fluidDiag(mgrid.A, mgrid.L,i,j,k,mgrid.dimensions,subcell);
			float e = diag - (left*left) - (bottom*bottom) - (back*back);
			if(diag>0){
				if( e < a*diag ){
					e = diag;
				}
				pc->setCell(i,j,k, 1.0f/sqrt(e));
			}
		}
	}
}

void flipDivergence(macgrid& mgrid){
	int x = (int)mgrid.dimensions.x; int y = (int)mgrid.dimensions.y; int z = (int)mgrid.dimensions.z;
	// #pragma omp parallel for
	for( int gn=0; gn<x*y*z; gn++ ) { 
		int i=(gn%((x)*(y)))%(z); int j=(gn%((x)*(y)))/(z); int k = gn/((x)*(y)); 
		float flippedD = -mgrid.D->getCell(i,j,k);
		mgrid.D->setCell(i,j,k,flippedD);
	}
}

void solve(macgrid& mgrid, const int& subcell){
	//flip divergence
	cout << "Flipping divergence..." << endl;
	flipDivergence(mgrid);

	//build preconditioner
	cout << "Building preconditioner matrix..." << endl;
	floatgrid* preconditioner = new floatgrid(0.0f);
	buildPreconditioner(preconditioner, mgrid, subcell);

	delete preconditioner;
}

}

#endif
