// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: geom.hpp
// Classes for geom frame operations and stuff

#ifndef GEOM_HPP
#define GEOM_HPP

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#define SHARED __shared__
#else
#define HOST
#define DEVICE
#define SHARED
#endif

#include "../spatial/bvh.hpp"
#include "../ray/ray.hpp"

enum GeomType{MESH, VOLUME, CURVE, NONE};

namespace geomCore {

//====================================
// Class Declarations
//====================================

class Geom;
class GeomFrame;
class GeomInterface;

class Geom {
	public:
		HOST DEVICE Geom();
		HOST DEVICE Geom(GeomInterface* geom);
		HOST DEVICE ~Geom();

		HOST DEVICE void SetContents(GeomInterface* geom);
        HOST DEVICE void Intersect(const rayCore::Ray& r, 
        						   spaceCore::TraverseAccumulator& result);
        HOST DEVICE GeomType GetType();

		GeomInterface*  	    m_geom;
		unsigned int    		m_id;
};

class GeomInterface {
	public:
		HOST DEVICE GeomInterface(){};
		HOST DEVICE ~GeomInterface(){};

		HOST DEVICE virtual void Intersect(const rayCore::Ray& r,
										   spaceCore::TraverseAccumulator& result) = 0;
		HOST DEVICE virtual GeomType GetType() = 0;
		HOST DEVICE virtual unsigned int GetID() = 0;
};

class GeomFrame {
	public:
		GeomFrame();
		GeomFrame(const glm::vec3& t, const glm::vec3& r, const glm::vec3& s);
		~GeomFrame();

		void SetContents(const glm::vec3& t, const glm::vec3& r, const glm::vec3& s);

		glm::vec3                   m_translation;
		glm::vec3                   m_rotation;
		glm::vec3                   m_scale;
        unsigned int                m_id;
};

}

#endif
