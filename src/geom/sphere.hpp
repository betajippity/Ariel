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
	
class Sphere: public Geom {
	public:
		//Initializers
		Sphere();
		Sphere(int subdivCount);
		~Sphere();

		//Getters
		objCore::objContainer* Tesselate();
		objCore::objContainer* Tesselate(const glm::vec3& center, const float& radius);

		//Data
		int m_subdivs;

	private:
		glm::vec3 GetPointOnSphereByAngles(float angle1, float angle2);
};
}

#endif
