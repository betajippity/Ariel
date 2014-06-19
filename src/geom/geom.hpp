// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: geom.hpp
// Defines the generic geometry pure virtual class that all geometry will inherit from

#ifndef GEOM_HPP
#define GEOM_HPP

#include <glm/glm.hpp>
#include "obj/obj.hpp"

namespace geomCore {
//====================================
// Class Declarations
//====================================

class Geom{
	public:
		//Initializers
		Geom(){};
		virtual ~Geom(){};

		//Interactions
		virtual objCore::Obj* Tesselate() = 0;
};
}

#endif