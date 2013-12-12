// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: sphere.hpp
// Sphere class, meant for fluid initialization and barrier placement. Implements geom.hpp

#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "geom.hpp"

using namespace std;
using namespace glm;

namespace geomCore {
//====================================
// Class Declarations
//====================================

class sphere: public geom {
	public:
		sphere(const float& radius, const vec3& center, const geomtype& type);
		sphere();
		~sphere();

		bool isPointInside(const vec3& point, const float& scale);
		bool isPointInsideWithThickness(const vec3& point, const float& thickness, const float& scale);
		geomtype getType();

	private:
		float radius;
		vec3 center;
		// geomtype type;
};
}

#endif
