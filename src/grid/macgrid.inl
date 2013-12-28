// Ariel: FLIP Fluid Simulator
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

enum geomtype {SOLID=2, FLUID=1, AIR=0};

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

	floatgrid* D; //divergence 
	floatgrid* P; //pressure
	intgrid* A; //cell type
	floatgrid* L; //internal lightweight SDF for project step

	gridtype type;
};

struct particle{
	vec3 p; //position
	vec3 u; //velocity
	vec3 n; //normal
	float density;
	float mass;
	int type;
	vec3 t;
	vec3 t2;
	bool invalid;
};

//Forward declarations for externed inlineable methods
extern inline particle createParticle(const vec3& position, const vec3& velocity, const vec3& normal, 
									  const float& density);
extern inline macgrid createMacgrid(const vec3& dimensions, const gridtype& type);
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

macgrid createMacgrid(const vec3& dimensions, const gridtype& type){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	macgrid m;
	m.dimensions = dimensions;
	m.u_x = new floatgrid(type, vec3(x+1,y,z), 0.0f);
	m.u_y = new floatgrid(type, vec3(x,y+1,z), 0.0f);
	m.u_z = new floatgrid(type, vec3(x,y,z+1), 0.0f);
	m.D = new floatgrid(type, vec3(x,y,z), 0.0f);
	m.P = new floatgrid(type, vec3(x,y,z), 0.0f);
	m.A = new intgrid(type, vec3(x,y,z), 0);
	m.L = new floatgrid(type, vec3(x,y,z), 1.6f);
	m.type = type;
	return m;
}

void clearMacgrid(macgrid& m){
	delete m.u_x;
	delete m.u_y;
	delete m.u_z;
	delete m.D;
	delete m.P;
}
}

#endif
