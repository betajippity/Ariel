// ObjCore4.0: An (improved) obj mesh wrangling library. Part of TAKUA Render.
// Written by Yining Karl Li
//
// File: obj.hpp
// A simple obj class. This library superscedes ObjCore 1.x/2.x. 

#ifndef OBJ_HPP
#define OBJ_HPP

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#else
#define HOST
#define DEVICE
#endif

#include <sstream>
#include <fstream>
#include <iostream>
#include "../../utilities/utilities.h"

namespace objCore {

//====================================
// Struct Declarations
//====================================

struct Point {
    glm::vec3 m_position;
    glm::vec3 m_normal;
    glm::vec2 m_uv;

    Point(){}

    Point(const glm::vec3& p, const glm::vec3& n, const glm::vec2& u){
    	m_position = p;
    	m_normal = n;
    	m_uv = u;
    }
};
    
struct Poly {
    Point m_vertex0;
    Point m_vertex1;
    Point m_vertex2;
    Point m_vertex3;
    unsigned int m_id;

    Poly(){}

    //Builds a quad struct with the given data
    Poly(const Point& v0, const Point& v1, const Point& v2, const Point& v3){
    	m_vertex0 = v0;
    	m_vertex1 = v1;
    	m_vertex2 = v2;
    	m_vertex3 = v3;
    }

    Poly(const Point& v0, const Point& v1, const Point& v2, const Point& v3, 
         const unsigned int& id){
    	m_vertex0 = v0;
    	m_vertex1 = v1;
    	m_vertex2 = v2;
    	m_vertex3 = v3;
        m_id = id;
    }

    //Builds a triangle struct with the given data
    Poly(const Point& v0, const Point& v1, const Point& v2){
    	m_vertex0 = v0;
    	m_vertex1 = v1;
    	m_vertex2 = v2;
    	m_vertex3 = v0; //in a triangle, the fourth vertex is the same as the first
    }

    Poly(const Point& v0, const Point& v1, const Point& v2, const unsigned int& id){
    	m_vertex0 = v0;
    	m_vertex1 = v1;
    	m_vertex2 = v2;
    	m_vertex3 = v0; //in a triangle, the fourth vertex is the same as the firsit
        m_id = id;
    }
};

//====================================
// Class Declarations
//====================================

class Obj {
	public:
		Obj();
		Obj(const std::string& filename);
		~Obj();

		void BakeTransform(const glm::mat4& transform);

		bool ReadObj(const std::string& filename);
		bool WriteObj(const std::string& filename);

		HOST DEVICE Poly GetPoly(const unsigned int& polyIndex);
		HOST DEVICE static Poly TransformPoly(const Poly& p, const glm::mat4& m);
		HOST DEVICE static Point TransformPoint(const Point& p, const glm::mat4& m);

		unsigned int    m_numberOfVertices;
	    glm::vec3*      m_vertices;
	    unsigned int    m_numberOfNormals;
	    glm::vec3*      m_normals;
	    unsigned int    m_numberOfUVs;
	    glm::vec2*      m_uvs;
	    unsigned int    m_numberOfPolys;
	    glm::uvec4*     m_polyVertexIndices;
	    glm::uvec4*     m_polyNormalIndices;
	    glm::uvec4*     m_polyUVIndices;
        unsigned int    m_id;
        bool            m_keep;
        
	private:
		void PrereadObj(const std::string& filename);

};
}

#endif
