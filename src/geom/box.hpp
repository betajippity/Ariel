// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: box.hpp
// Box class, meant for fluid initialization and barrier placement. Implements geom.hpp

#ifndef BOX_HPP
#define BOX_HPP

#include "geom.hpp"

using namespace std;
using namespace glm;

namespace geomCore {
//====================================
// Class Declarations
//====================================

class box {
	public:
		box(const vec3& lowerCorner, const vec3& upperCorner, const geomtype& type);
		box();
		~box();

		bool isParticleInside(const vec3& point);
		bool isParticleInsideWithThickness(const vec3& point, const float& thickness);
		geomtype getType();

	private:
		vec3 upperCorner;
		vec3 lowerCorner;
		geomtype type;
};
}

#endif
