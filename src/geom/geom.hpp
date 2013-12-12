// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: geom.hpp
// Defines the generic geometry pure virtual class that all geometry will inherit from

#ifndef GEOM_HPP
#define GEOM_HPP

#include <glm/glm.hpp>
#include "obj/objcontainer.hpp"

using namespace std;
using namespace glm;

namespace geomCore {
//====================================
// Class Declarations
//====================================

class geom{
	public:
		//Initializers
		geom(){};
		virtual ~geom(){};

		//Interactions
		virtual objCore::objContainer* tesselate() = 0;
};
}

#endif