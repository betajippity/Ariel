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

class box: public geom {
	public:
		box(const vec3& lowerCorner, const vec3& upperCorner, const geomtype& type);
		box();
		~box();

		bool isPointInside(const vec3& point, const float& scale);
		bool isPointInsideWithThickness(const vec3& point, const float& thickness, const float& scale);
		geomtype getType();

	private:
		vec3 upperCorner;
		vec3 lowerCorner;
		geomtype type;
};
}

#endif
