// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: solver.inl
// Breakout file for PCG solver

#ifndef SOLVER_INL
#define SOLVER_INL

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
extern inline void Solve(MacGrid& mgrid, const int& subcell, const bool& verbose);
inline void BuildPreconditioner(Grid<float>* pc, MacGrid& mgrid, int subcell);
inline void SolveConjugateGradient(MacGrid& mgrid, Grid<float>* pc, int subcell, 
								   const bool& verbose);
inline void ComputeAx(Grid<int>* A, Grid<float>* L, Grid<float>* X, Grid<float>* target, 
					  glm::vec3 dimensions, int subcell);
inline float XRef(Grid<int>* A, Grid<float>* L, Grid<float>* X, glm::vec3 f, glm::vec3 p, 
				  glm::vec3 dimensions, int subcell);
inline void Op(Grid<int>* A, Grid<float>* X, Grid<float>* Y, Grid<float>* target, float alpha, 
			   glm::vec3 dimensions);
inline float Product(Grid<int>* A, Grid<float>* X, Grid<float>* Y, glm::vec3 dimensions);
inline void ApplyPreconditioner(Grid<float>* Z, Grid<float>* R, Grid<float>* P, Grid<float>* L, 
								Grid<int>* A, glm::vec3 dimensions);

//====================================
// Function Implementations
//====================================

//Takes a grid, multiplies everything by -1
void FlipGrid(Grid<float>* grid, glm::vec3 dimensions){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){
		  		for(unsigned int j=0; j<y; ++j){
					for(unsigned int k=0; k<z; ++k){
						float flipped = -grid->GetCell(i,j,k);
						grid->SetCell(i,j,k,flipped);
					}
				}      	
	      	}
	    }
    );
}

//Helper for preconditioner builder
float ARef(Grid<int>* A, int i, int j, int k, int qi, int qj, int qk, glm::vec3 dimensions){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	if( i<0 || i>x-1 || j<0 || j>y-1 || k<0 || k>z-1 || A->GetCell(i,j,k)!=FLUID ){ //if not liquid
		return 0.0;
	} 
	//if not liquid
	if( qi<0 || qi>x-1 || qj<0 || qj>y-1 || qk<0 || qk>z-1 || A->GetCell(qi,qj,qk)!=FLUID ){ 
		return 0.0;
	} 
	return -1.0;	
}

//Helper for preconditioner builder
float PRef(Grid<float>* p, int i, int j, int k, glm::vec3 dimensions){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	if( i<0 || i>x-1 || j<0 || j>y-1 || k<0 || k>z-1 || p->GetCell(i,j,k)!=FLUID ){ //if not liquid
		return 0.0f;
	} 
	return p->GetCell(i,j,k);
}

//Helper for preconditioner builder
float ADiag(Grid<int>* A, Grid<float>* L, int i, int j, int k, glm::vec3 dimensions, int subcell){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	float diag = 6.0;
	if( A->GetCell(i,j,k) != FLUID ){
		return diag;
	}
	int q[][3] = { {i-1,j,k}, {i+1,j,k}, {i,j-1,k}, {i,j+1,k}, {i,j,k-1}, {i,j,k+1} };
	for( int m=0; m<6; m++ ) {
		int qi = q[m][0];
		int qj = q[m][1];
		int qk = q[m][2];
		if( qi<0 || qi>x-1 || qj<0 || qj>y-1 || qk<0 || qk>z-1 ){
			diag -= 1.0;
		}
		else if( A->GetCell(qi,qj,qk)==AIR && subcell ) {
			diag -= L->GetCell(qi,qj,qk)/glm::min(1.0e-6f,L->GetCell(i,j,k));
		}
	}
	
	return diag;
}

//Does what it says
void BuildPreconditioner(Grid<float>* pc, MacGrid& mgrid, int subcell){
	int x = (int)mgrid.m_dimensions.x; int y = (int)mgrid.m_dimensions.y; 
	int z = (int)mgrid.m_dimensions.z;
	float a = 0.25f;
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){
				for(unsigned int j=0; j<y; ++j){
					for(unsigned int k=0; k<z; ++k){
						if(mgrid.m_A->GetCell(i,j,k)==FLUID){	
							float left = ARef(mgrid.m_A,i-1,j,k,i,j,k,mgrid.m_dimensions) * 
										 PRef(pc,i-1,j,k,mgrid.m_dimensions);
							float bottom = ARef(mgrid.m_A,i,j-1,k,i,j,k,mgrid.m_dimensions) * 
										   PRef(pc,i,j-1,k,mgrid.m_dimensions);
							float back = ARef(mgrid.m_A,i,j,k-1,i,j,k,mgrid.m_dimensions) * 
										 PRef(pc,i,j,k-1,mgrid.m_dimensions);
							float diag = ADiag(mgrid.m_A, mgrid.m_L,i,j,k,mgrid.m_dimensions,
											   subcell);
							float e = diag - (left*left) - (bottom*bottom) - (back*back);
							if(diag>0){
								if( e < a*diag ){
									e = diag;
								}
								pc->SetCell(i,j,k, 1.0f/glm::sqrt(e));
							}
						}
					}
				}
			}
		}
	);
}

//Helper for PCG solver: read X with clamped bounds
float XRef(Grid<int>* A, Grid<float>* L, Grid<float>* X, glm::vec3 f, glm::vec3 p, 
		   glm::vec3 dimensions, int subcell){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	int i = glm::min(glm::max(0,(int)p.x),x-1); int fi = (int)f.x;
	int j = glm::min(glm::max(0,(int)p.y),y-1); int fj = (int)f.y;
	int k = glm::min(glm::max(0,(int)p.z),z-1); int fk = (int)f.z;
	if(A->GetCell(i,j,k) == FLUID){
		return X->GetCell(i,j,k);
	}else if(A->GetCell(i,j,k) == SOLID){
		return X->GetCell(fi,fj,fk);
	} 
	if(subcell){
		return L->GetCell(i,j,k)/glm::min(1.0e-6f,L->GetCell(fi,fj,fk))*X->GetCell(fi,fj,fk);
	}else{
		return 0.0f;
	}
}

// target = X + alpha*Y
void Op(Grid<int>* A, Grid<float>* X, Grid<float>* Y, Grid<float>* target, float alpha, 
		glm::vec3 dimensions){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	for(unsigned int j=0; j<y; ++j){
		for(unsigned int k=0; k<z; ++k){
			//this parallel loop has to be the inner loop or else MSVC will barf
			tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
				[=](const tbb::blocked_range<unsigned int>& r){
					for(unsigned int i=r.begin(); i!=r.end(); ++i){
						if(A->GetCell(i,j,k)==FLUID){
							float targetval = X->GetCell(i,j,k)+alpha*Y->GetCell(i,j,k);
							target->SetCell(i,j,k,targetval);
						}else{
							target->SetCell(i,j,k,0.0f);
						}				
					}
				}
			);
		}
	}
}

// ans = x^T * x
float Product(Grid<int>* A, Grid<float>* X, Grid<float>* Y, glm::vec3 dimensions){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	float result = 0.0f;
	for(unsigned int i=0; i<x; i++){
		for(unsigned int j=0; j<y; j++){
			for(unsigned int k=0; k<z; k++){
				if(A->GetCell(i,j,k)==FLUID){
					result += X->GetCell(i,j,k) * Y->GetCell(i,j,k);
				}
			}
		}
	}
	return result;
}

//Helper for PCG solver: target = AX
void ComputeAx(Grid<int>* A, Grid<float>* L, Grid<float>* X, Grid<float>* target, 
			   glm::vec3 dimensions, int subcell){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	float n = (float)glm::max(glm::max(x,y),z);
	float h = 1.0f/(n*n);
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){
				for(unsigned int j=0; j<y; ++j){
					for(unsigned int k=0; k<z; ++k){
						if(A->GetCell(i,j,k) == FLUID){
							float result = (6.0f*X->GetCell(i,j,k)
											-XRef(A, L, X, glm::vec3(i,j,k), glm::vec3(i+1,j,k), 
												  dimensions, subcell)
											-XRef(A, L, X, glm::vec3(i,j,k), glm::vec3(i-1,j,k), 
												  dimensions, subcell)
											-XRef(A, L, X, glm::vec3(i,j,k), glm::vec3(i,j+1,k), 
												  dimensions, subcell)
											-XRef(A, L, X, glm::vec3(i,j,k), glm::vec3(i,j-1,k), 
												  dimensions, subcell)
											-XRef(A, L, X, glm::vec3(i,j,k), glm::vec3(i,j,k+1), 
												  dimensions, subcell)
											-XRef(A, L, X, glm::vec3(i,j,k), glm::vec3(i,j,k-1), 
												  dimensions, subcell)
											)/h;
							target->SetCell(i,j,k,result);
						} else {
							target->SetCell(i,j,k,0.0f);
						}
					}
				}
			}
		}
	);
}

void ApplyPreconditioner(Grid<float>* Z, Grid<float>* R, Grid<float>* P, Grid<float>* L, 
						 Grid<int>* A, glm::vec3 dimensions){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	Grid<float>* Q = new Grid<float>(dimensions, 0.0f);

	// LQ = R
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){
				for(unsigned int j=0; j<y; ++j){
					for(unsigned int k=0; k<z; ++k){
						if(A->GetCell(i,j,k) == FLUID) {
							float left = ARef(A,i-1,j,k,i,j,k,dimensions)*
										 PRef(P,i-1,j,k,dimensions)*PRef(Q,i-1,j,k,dimensions);
							float bottom = ARef(A,i,j-1,k,i,j,k,dimensions)*
										   PRef(P,i,j-1,k,dimensions)*PRef(Q,i,j-1,k,dimensions);
							float back = ARef(A,i,j,k-1,i,j,k,dimensions)*
										 PRef(P,i,j,k-1,dimensions)*PRef(Q,i,j,k-1,dimensions);
							
							float t = R->GetCell(i,j,k) - left - bottom - back;
							float qVal = t * P->GetCell(i,j,k);
							Q->SetCell(i,j,k,qVal);
						}
					}
				}
			}
		}
	);

	// L^T Z = Q
	for(int j=y-1; j>=0; j--){
		for(int k=z-1; k>=0; k--){
			//this parallel loop has to be the inner loop or else MSVC will barf
			tbb::parallel_for(tbb::blocked_range<int>(-1,x-1),
				[=](const tbb::blocked_range<int>& r){
					for(int i=r.end(); i!=r.begin(); i--){
						if(A->GetCell(i,j,k) == FLUID){
							float right = ARef(A,i,j,k,i+1,j,k,dimensions)*
										  PRef(P,i,j,k,dimensions)*PRef(Z,i+1,j,k,dimensions);
							float top = ARef(A,i,j,k,i,j+1,k,dimensions)*
										PRef(P,i,j,k,dimensions)*PRef(Z,i,j+1,k,dimensions);
							float front = ARef(A,i,j,k,i,j,k+1,dimensions)*
										  PRef(P,i,j,k,dimensions)*PRef(Z,i,j,k+1,dimensions);
						
							float t = Q->GetCell(i,j,k) - right - top - front;
							float zVal = t * P->GetCell(i,j,k);
							Z->SetCell(i,j,k,zVal);
						}
					}
				}
			);
		}
	}
	delete Q;
}

//Does what it says
void SolveConjugateGradient(MacGrid& mgrid, Grid<float>* PC, int subcell, const bool& verbose){
	int x = (int)mgrid.m_dimensions.x; int y = (int)mgrid.m_dimensions.y; 
	int z = (int)mgrid.m_dimensions.z;

	Grid<float>* R = new Grid<float>(mgrid.m_dimensions, 0.0f);
	Grid<float>* Z = new Grid<float>(mgrid.m_dimensions, 0.0f);
	Grid<float>* S = new Grid<float>(mgrid.m_dimensions, 0.0f);

	//note: we're calling pressure "mgrid.P" instead of x

	ComputeAx(mgrid.m_A, mgrid.m_L, mgrid.m_P, Z, mgrid.m_dimensions, subcell);	// z = apply A(x)
	Op(mgrid.m_A, mgrid.m_D, Z, R, -1.0f, mgrid.m_dimensions);                // r = b-Ax
	float error0 = Product(mgrid.m_A, R, R, mgrid.m_dimensions);			// error0 = product(r,r)

	// z = f(r), aka preconditioner step
	ApplyPreconditioner(Z, R, PC, mgrid.m_L, mgrid.m_A, mgrid.m_dimensions);	

	//s = z. TODO: replace with VDB deep copy?

	for(unsigned int j=0; j<y; ++j ){
		for(unsigned int k=0; k<z; ++k ){
			//this parallel loop has to be the inner loop or else MSVC will barf
			tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
				[=](const tbb::blocked_range<unsigned int>& r){
					for(unsigned int i=r.begin(); i!=r.end(); ++i){
						S->SetCell(i,j,k,Z->GetCell(i,j,k));
					}
				}
			);	
		}
	}


	float eps = 1.0e-2f * (x*y*z);
	float a = Product(mgrid.m_A, Z, R, mgrid.m_dimensions);					// a = product(z,r)

	for( int k=0; k<x*y*z; k++){
		//Solve current iteration
		ComputeAx(mgrid.m_A, mgrid.m_L, S, Z, mgrid.m_dimensions, subcell);	// z = applyA(s)
		float alpha = a/Product(mgrid.m_A, Z, S, mgrid.m_dimensions);		// alpha = a/(z . s)
		Op(mgrid.m_A, mgrid.m_P, S, mgrid.m_P, alpha, mgrid.m_dimensions);		// x = x + alpha*s
		Op(mgrid.m_A, R, Z, R, -alpha, mgrid.m_dimensions);					// r = r - alpha*z;
		float error1 = Product(mgrid.m_A, R, R, mgrid.m_dimensions);		// error1 = product(r,r)
        error0 = glm::max(error0, error1);
        //Output progress
        float rate = 1.0f - glm::max(0.0f,glm::min(1.0f,(error1-eps)/(error0-eps)));
        if(verbose){
        	std::cout << "PCG Iteration " << k+1 << ": " << 100.0f*pow(rate,6) << "%% solved" 
        			  << std::endl;
        }
        if(error1<=eps){
        	break;
        }
        //Prep next iteration
        // z = f(r)
        ApplyPreconditioner(Z, R, PC, mgrid.m_L, mgrid.m_A, mgrid.m_dimensions);
        float a2 = Product(mgrid.m_A, Z, R, mgrid.m_dimensions);				// a2 = product(z,r)
        float beta = a2/a;													// beta = a2/a
        Op(mgrid.m_A, Z, S, S, beta, mgrid.m_dimensions);						// s = z + beta*s
        a = a2;
	}

	delete R;
	delete Z;
	delete S;
}

void Solve(MacGrid& mgrid, const int& subcell, const bool& verbose){

	//if in VDB mode, force to single threaded to prevent VDB write issues. 
	//this is a kludgey fix for now.
	// if(mgrid.type==VDB){
	// 	omp_set_num_threads(1);
	// }

	//flip divergence
	FlipGrid(mgrid.m_D, mgrid.m_dimensions);

	//build preconditioner
	Grid<float>* preconditioner = new Grid<float>(mgrid.m_dimensions, 0.0f);
	BuildPreconditioner(preconditioner, mgrid, subcell);

	//solve conjugate gradient
	SolveConjugateGradient(mgrid, preconditioner, subcell, verbose);

	delete preconditioner;

	// if(mgrid.type==VDB){
	// 	omp_set_num_threads(omp_get_num_procs());
	// }
}
}

#endif
