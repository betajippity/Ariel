// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: sphere.hpp
// Defines the unit sphere geometry class, inherits from the generic geom class 

#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "geom.hpp"

namespace geomCore {
//====================================
// Class Declarations
//====================================
	
class sphere: public geom {
	public:
		//Initializers
		sphere();
		sphere(int subdivCount);
		~sphere();

		//Getters
		objCore::objContainer* tesselate();
		objCore::objContainer* tesselate(const glm::vec3& center, const float& radius);

		//Data
		int subdivs;

	private:
		glm::vec3 getPointOnSphereByAngles(float angle1, float angle2);
};
}

#endif
