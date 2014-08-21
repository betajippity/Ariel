// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: macgrid.inl
// MAC grid built on top of a bunch of floatgrids

#ifndef MACGRID_INL
#define MACGRID_INL

#include "grid.hpp"
#include "../utilities/utilities.h"

enum geomtype {SOLID=2, FLUID=1, AIR=0};

namespace fluidCore {
//====================================
// Struct and Function Declarations
//====================================

//A macgrid in this simulator is built up entirely out of VDB volumes
struct MacGrid{
	glm::vec3		m_dimensions;

	//face velocities
	Grid<float>*	m_u_x;
	Grid<float>*	m_u_y;
	Grid<float>*	m_u_z; 
	//technically this is the part that is an actual MAC grid, the rest is other useful stuff

	Grid<float>*	m_D; //divergence 
	Grid<float>*	m_P; //pressure
	Grid<int>*		m_A; //cell type
	Grid<float>*	m_L; //internal lightweight SDF for project step
};

struct Particle{
	glm::vec3		m_p; //position
	glm::vec3		m_u; //velocity
	glm::vec3		m_n; //normal
	float			m_density;
	float			m_mass;
	int				m_type;
	glm::vec3		m_t;
	glm::vec3		m_t2;
    glm::vec3       m_pt; //copy of previous position used for bound checks
	bool			m_invalid;
	bool            m_temp2;
    bool			m_temp;
};

//Forward declarations for externed inlineable methods
extern inline Particle CreateParticle(const glm::vec3& position, const glm::vec3& velocity, 
									  const glm::vec3& normal, const float& density);
extern inline MacGrid CreateMacgrid(const glm::vec3& dimensions);
extern inline void ClearMacgrid(MacGrid& m);

//====================================
// Function Implementations
//====================================

Particle CreateParticle(const glm::vec3& position, const glm::vec3& velocity, 
						const glm::vec3& normal, const float& density){
	Particle p;
	p.m_p = position;
	p.m_u = velocity;
	p.m_n = normal;
	p.m_density= density;
	return p;
}

MacGrid CreateMacgrid(const glm::vec3& dimensions){
	int x = (int)dimensions.x; int y = (int)dimensions.y; int z = (int)dimensions.z;
	MacGrid m;
	m.m_dimensions = dimensions;
	m.m_u_x = new Grid<float>(glm::vec3(x+1,y,z), 0.0f);
	m.m_u_y = new Grid<float>(glm::vec3(x,y+1,z), 0.0f);
	m.m_u_z = new Grid<float>(glm::vec3(x,y,z+1), 0.0f);
	m.m_D = new Grid<float>(glm::vec3(x,y,z), 0.0f);
	m.m_P = new Grid<float>(glm::vec3(x,y,z), 0.0f);
	m.m_A = new Grid<int>(glm::vec3(x,y,z), 0);
	m.m_L = new Grid<float>(glm::vec3(x,y,z), 1.6f);
	return m;
}

void ClearMacgrid(MacGrid& m){
	delete m.m_u_x;
	delete m.m_u_y;
	delete m.m_u_z;
	delete m.m_D;
	delete m.m_P;
}
}

#endif
