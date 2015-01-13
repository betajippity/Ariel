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
#include "../../ray/ray.hpp"
#include "../../spatial/aabb.hpp"

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

        //Standard interface functions for accel. structures
        HOST DEVICE rayCore::Intersection IntersectElement(const unsigned int& primID,
                                                           const rayCore::Ray& r);
        HOST DEVICE spaceCore::Aabb GetElementAabb(const unsigned int& primID);
        HOST DEVICE unsigned int GetNumberOfElements();

        HOST DEVICE static inline rayCore::Intersection RayTriangleTest(const glm::vec3& v0,
                                                                        const glm::vec3& v1,
                                                                        const glm::vec3& v2,
                                                                        const glm::vec3& n0,
                                                                        const glm::vec3& n1,
                                                                        const glm::vec3& n2,
                                                                        const glm::vec2& u0,
                                                                        const glm::vec2& u1,
                                                                        const glm::vec2& u2,
                                                                        const rayCore::Ray& r);

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
        HOST DEVICE rayCore::Intersection TriangleTest(const unsigned int& polyIndex, 
                                                       const rayCore::Ray& r, 
                                                       const bool& checkQuad);
        HOST DEVICE rayCore::Intersection QuadTest(const unsigned int& polyIndex,
                                                   const rayCore::Ray& r);
};

class InterpolatedObj {
    public:
        InterpolatedObj();
        InterpolatedObj(objCore::Obj* obj0, objCore::Obj* obj1);
        ~InterpolatedObj();

        HOST DEVICE Poly GetPoly(const unsigned int& polyIndex, const float& interpolation);

        //Standard interface functions for accel. structures
        HOST DEVICE rayCore::Intersection IntersectElement(const unsigned int& primID,
                                                           const rayCore::Ray& r);
        HOST DEVICE spaceCore::Aabb GetElementAabb(const unsigned int& primID);
        HOST DEVICE unsigned int GetNumberOfElements();

        objCore::Obj*   m_obj0;
        objCore::Obj*   m_obj1;

    private:
        HOST DEVICE rayCore::Intersection TriangleTest(const unsigned int& polyIndex, 
                                                       const rayCore::Ray& r, 
                                                       const bool& checkQuad);
        HOST DEVICE rayCore::Intersection QuadTest(const unsigned int& polyIndex,
                                                   const rayCore::Ray& r); 
};
}

#endif
