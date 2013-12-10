// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: geom.hpp
// Geometry pure virtual base class, meant for fluid initialization and barrier placement

#ifndef GEOM_HPP
#define GEOM_HPP

#include "../utilities/utilities.h"

using namespace std;
using namespace glm;

enum geomtype {SOLID, FLUID};

namespace geomCore {
//====================================
// Class Declarations
//====================================

class geom {
	public:
		geom(){};
		virtual ~geom(){};

		virtual bool isParticleInside(const vec3& point) = 0;
		virtual bool isParticleInsideWithThickness(const vec3& point, const float& thickness) = 0;
		virtual geomtype getType() = 0;
};
}

#endif
