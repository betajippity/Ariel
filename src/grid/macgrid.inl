// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: macgrid.inl
// MAC grid built on top of a bunch of floatgrids

#ifndef MACGRID_INL
#define MACGRID_INL

#include "floatgrid.hpp"
#include "intgrid.hpp"
#include "../utilities/utilities.h"

using namespace std;
using namespace glm;

#define AIR 0
#define FLUID 1
#define SOLID 2

namespace fluidCore {
//====================================
// Struct and Function Declarations
//====================================

//A macgrid in this simulator is built up entirely out of VDB volumes
struct macgrid{
	vec3 dimensions;

	//face velocities
	floatgrid* u_x;
	floatgrid* u_y;
	floatgrid* u_z; //technically this is the part that is an actual MAC grid, the rest is other useful stuff

	intgrid* A; //stores cell type: air, fluid, or solid
	floatgrid* L; //level set for project step
	floatgrid* P; //pressure
};

struct particle{
	vec3 p; //position
	vec3 u; //velocity
	vec3 n; //normal
	float density;
	float mass;
	int type;
};

//Forward declarations for externed inlineable methods
extern inline particle createParticle(const vec3& position, const vec3& velocity, const vec3& normal, 
									  const float& density);
extern inline macgrid createMacgrid(const vec3& dimensions);
extern inline void clearMacgrid(macgrid& m);

//====================================
// Function Implementations
//====================================

particle createParticle(const vec3& position, const vec3& velocity, const vec3& normal, 
						const float& density){
	particle p;
	p.p = position;
	p.u = velocity;
	p.n = normal;
	p.density= density;
	return p;
}

macgrid createMacgrid(const vec3& dimensions){
	macgrid m;
	m.dimensions = dimensions;
	m.u_x = new floatgrid(0.0f);
	m.u_y = new floatgrid(0.0f);
	m.u_z = new floatgrid(0.0f);
	m.A = new intgrid(0);
	m.L = new floatgrid(0.0f);
	m.P = new floatgrid(0.0f);
	return m;
}

void clearMacgrid(macgrid& m){
	delete m.u_x;
	delete m.u_y;
	delete m.u_z;
	delete m.A;
	delete m.L;
	delete m.P;
}
}

#endif
