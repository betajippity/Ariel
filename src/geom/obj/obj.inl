// ObjCore2.5: An (improved) obj mesh wrangling library. Part of TAKUA Render.
// Written by Yining Karl Li
//
// File: obj.inl
// Implements data struct and inlineable functions for ObjCore2.5

#ifndef OBJ_INL
#define OBJ_INL

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#else
#define HOST
#define DEVICE
#endif

#include <glm/glm.hpp>
#include "../../utilities/utilities.h"
 
using namespace std;
using namespace glm;

namespace objCore {
//====================================
// Struct and Function Declarations
//====================================
    
/*Defines all data that make up an obj. Only stores triangles and quads*/
struct obj {
    int numberOfVertices;
    vec3* vertices;
    int numberOfNormals;
    vec3* normals;
    int numberOfUVs;
    vec2* uvs;
    int numberOfPolys;
    vec4* polyVertexIndices;
    vec4* polyNormalIndices;
    vec4* polyUVIndices;
};
  
struct point {
    vec3 position;
    vec3 normal;
    vec2 uv;
};
    
struct poly {
    point vertex0;
    point vertex1;
    point vertex2;
    point vertex3;
};
    
//Forward declarations for externed inlineable methods
extern inline obj* createObj(int numberOfVertices, vec3* vertices, int numberOfNormals, vec3* normals, 
                             int numberOfUVs, vec2* uvs, int numberOfPolys, vec4* polyVertexIndices, 
                             vec4* polyNormalIndices, vec4* polyUVIndices);
extern inline void clearObj(obj* mesh);
HOST DEVICE extern inline point createPoint(vec3 position, vec3 normal, vec2 uv);
HOST DEVICE extern inline poly createPoly(point v1, point v2, point v3);
HOST DEVICE extern inline poly createPoly(point v1, point v2, point v3, point v4);
HOST DEVICE extern inline poly transformPoly(poly p, mat4 m);
HOST DEVICE extern inline point transformPoint(point p, mat4 m);
HOST DEVICE extern inline poly getPoly(obj* mesh, int polyIndex);

//====================================
// Function Implementations
//====================================
    
/*Build an obj given inputs. No sanity-checking is done, this method trusts whatever it is given, so be
careful!*/
obj* createObj(int numberOfVertices, vec3* vertices, int numberOfNormals, vec3* normals, int numberOfUVs, 
               vec2* uvs, int numberOfPolys, vec4* polyVertexIndices, vec4* polyNormalIndices, 
               vec4* polyUVIndices){
    obj* mesh = new obj;
    mesh->numberOfVertices = numberOfVertices;
    mesh->vertices = vertices;
    mesh->numberOfNormals = numberOfNormals;
    mesh->normals = normals;
    mesh->numberOfUVs = numberOfUVs;
    mesh->uvs = uvs;
    mesh->numberOfPolys = numberOfPolys;
    mesh->polyVertexIndices = polyVertexIndices;
    mesh->polyNormalIndices = polyNormalIndices;
    mesh->polyUVIndices = polyUVIndices;
    return mesh;
}

//Deletes the contents of the obj struct
void clearObj(obj* mesh){
    delete [] mesh->vertices;
    delete [] mesh->normals;
    delete [] mesh->uvs;
    delete [] mesh->polyVertexIndices;
    delete [] mesh->polyNormalIndices;
    delete [] mesh->polyUVIndices;
}

//Builds a point struct with the given data
HOST DEVICE point createPoint(vec3 position, vec3 normal, vec2 uv){
    point p;
    p.position = position;
    p.normal = normalize(normal);
    p.uv = uv;
    return p;
}

//Builds a triangle struct with the given data
HOST DEVICE poly createPoly(point v1, point v2, point v3){
    poly t;
    t.vertex0 = v1;
    t.vertex1 = v2;
    t.vertex2 = v3;
    t.vertex3 = v1; //in a triangle, the fourth vertex is the same as the first
    return t;
}

//Builds a quad struct with the given data
HOST DEVICE poly createPoly(point v1, point v2, point v3, point v4){
    poly t;
    t.vertex0 = v1;
    t.vertex1 = v2;
    t.vertex2 = v3;
    t.vertex3 = v4;
    return t;
}

/*Return the requested face from the mesh, unless the index is out of range, in which case return a 
face of area zero*/
HOST DEVICE poly getPoly(obj* mesh, int polyIndex){
    point pNull = createPoint(vec3(0,0,0), vec3(0,1,0), vec2(0,0));
    if(polyIndex<0 || polyIndex>=mesh->numberOfPolys){
        return createPoly(pNull, pNull, pNull);
    }else{
        vec4 vertexIndices = mesh->polyVertexIndices[polyIndex];
        vec4 normalIndices = mesh->polyNormalIndices[polyIndex];
        vec4 uvIndices = mesh->polyUVIndices[polyIndex];
        point p1 = createPoint(mesh->vertices[(int)vertexIndices.x-1], 
                                       mesh->normals[(int)normalIndices.x-1], 
                                       mesh->uvs[(int)uvIndices.x-1]);
        point p2 = createPoint(mesh->vertices[(int)vertexIndices.y-1], 
                                       mesh->normals[(int)normalIndices.y-1], 
                                       mesh->uvs[(int)uvIndices.y-1]);
        point p3 = createPoint(mesh->vertices[(int)vertexIndices.z-1], 
                                       mesh->normals[(int)normalIndices.z-1], 
                                       mesh->uvs[(int)uvIndices.z-1]);
        point p4 = createPoint(mesh->vertices[(int)vertexIndices.w-1], 
                                       mesh->normals[(int)normalIndices.w-1], 
                                       mesh->uvs[(int)uvIndices.w-1]);
        if(vertexIndices.w-1<0){
            p4 = p1;
        }
        return createPoly(p1, p2, p3, p4);
    }
}

//Applies given transformation matrix to the given point
HOST DEVICE point transformPoint(point p, mat4 m){
    point r = p;
    r.position = vec3( utilityCore::multiply(m,vec4(p.position,1.0)) );
    r.normal = normalize(vec3( utilityCore::multiply(m,vec4(p.normal,0.0)) ));
    return r;
}

//Applies given transformation matrix to the given poly
HOST DEVICE poly transformPoly(poly p, mat4 m){
    poly r;
    r.vertex0 = transformPoint(p.vertex0, m);
    r.vertex1 = transformPoint(p.vertex1, m);
    r.vertex2 = transformPoint(p.vertex2, m);
    r.vertex3 = transformPoint(p.vertex3, m);
    return r;
}
}

#endif
