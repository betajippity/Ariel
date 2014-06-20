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
	
class Cube: public Geom {
	public:
		//Initializers
		Cube();
		~Cube();

		//Getters
		objCore::Obj* Tesselate();
		objCore::Obj* Tesselate(const glm::vec3& lowerCorner, const glm::vec3& upperCorner);
};
}

#endif