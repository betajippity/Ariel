// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: cube.hpp
// Defines the unit cube geometry class, inherits from the generic geom class 

#ifndef CUBE_HPP
#define CUBE_HPP

#include "geom.hpp"

using namespace std;
using namespace glm;

namespace geomCore {
//====================================
// Class Declarations
//====================================
	
class cube: public geom {
	public:
		//Initializers
		cube();
		~cube();

		//Getters
		objCore::objContainer* tesselate();
		objCore::objContainer* tesselate(const vec3& lowerCorner, const vec3& upperCorner);
};
}

#endif