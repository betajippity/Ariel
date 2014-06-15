// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: cube.hpp
// Defines the unit cube geometry class, inherits from the generic geom class 

#ifndef CUBE_HPP
#define CUBE_HPP

#include "geom.hpp"

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
    objCore::objContainer* tesselate(const glm::vec3& lowerCorner, const glm::vec3& upperCorner);
};
}

#endif