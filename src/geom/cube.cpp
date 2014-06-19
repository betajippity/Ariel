// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: cube.cpp
// Implements cube.hpp

#include <tbb/tbb.h>
#include "cube.hpp"
#include "../utilities/utilities.h" 

namespace geomCore{

//Default empty constructor defaults to 1 subdivs in axis and height
Cube::Cube(){	
}

//Boring empty destructor is boring
Cube::~Cube(){
}

//Literally just embed an obj because that's how we roll. Stupid but it works.
objCore::Obj* Cube::Tesselate(){
	glm::vec3* vertices = new glm::vec3[8];
	vertices[0] = glm::vec3(-0.5, -0.5,  0.5);
	vertices[1] = glm::vec3( 0.5, -0.5,  0.5);
	vertices[2] = glm::vec3(-0.5,  0.5,  0.5);
	vertices[3] = glm::vec3( 0.5,  0.5,  0.5);
	vertices[4] = glm::vec3(-0.5,  0.5, -0.5);
	vertices[5] = glm::vec3( 0.5,  0.5, -0.5);
	vertices[6] = glm::vec3(-0.5, -0.5, -0.5);
	vertices[7] = glm::vec3( 0.5, -0.5, -0.5);
	glm::vec3* normals = new glm::vec3[24];
	normals[ 0] = glm::vec3( 0.000000,  0.000000,  1.000000);
	normals[ 1] = glm::vec3( 0.000000,  0.000000,  1.000000);
	normals[ 2] = glm::vec3( 0.000000,  0.000000,  1.000000);
	normals[ 3] = glm::vec3( 0.000000,  0.000000,  1.000000);
	normals[ 4] = glm::vec3( 0.000000,  1.000000,  0.000000);
	normals[ 5] = glm::vec3( 0.000000,  1.000000,  0.000000);
	normals[ 6] = glm::vec3( 0.000000,  1.000000,  0.000000);
	normals[ 7] = glm::vec3( 0.000000,  1.000000,  0.000000);
	normals[ 8] = glm::vec3( 0.000000,  0.000000, -1.000000);
	normals[ 9] = glm::vec3( 0.000000,  0.000000, -1.000000);
	normals[10] = glm::vec3( 0.000000,  0.000000, -1.000000);
	normals[11] = glm::vec3( 0.000000,  0.000000, -1.000000);
	normals[12] = glm::vec3( 0.000000, -1.000000,  0.000000);
	normals[13] = glm::vec3( 0.000000, -1.000000,  0.000000);
	normals[14] = glm::vec3( 0.000000, -1.000000,  0.000000);
	normals[15] = glm::vec3( 0.000000, -1.000000,  0.000000);
	normals[16] = glm::vec3( 1.000000,  0.000000,  0.000000);
	normals[17] = glm::vec3( 1.000000,  0.000000,  0.000000);
	normals[18] = glm::vec3( 1.000000,  0.000000,  0.000000);
	normals[19] = glm::vec3( 1.000000,  0.000000,  0.000000);
	normals[20] = glm::vec3(-1.000000,  0.000000,  0.000000);
	normals[21] = glm::vec3(-1.000000,  0.000000,  0.000000);
	normals[22] = glm::vec3(-1.000000,  0.000000,  0.000000);
	normals[23] = glm::vec3(-1.000000,  0.000000,  0.000000);
	glm::vec2* uvs = new glm::vec2[4];
	uvs[0] = glm::vec2( 0, 0);
	uvs[1] = glm::vec2( 1, 0);
	uvs[2] = glm::vec2( 1, 1);
	uvs[3] = glm::vec2( 0, 1);
	glm::uvec4* polyVertexIndices = new glm::uvec4[6];
	polyVertexIndices[0] = glm::uvec4(1,2,4,3);
	polyVertexIndices[1] = glm::uvec4(3,4,6,5);
	polyVertexIndices[2] = glm::uvec4(5,6,8,7);
	polyVertexIndices[3] = glm::uvec4(7,8,2,1);
	polyVertexIndices[4] = glm::uvec4(2,8,6,4);
	polyVertexIndices[5] = glm::uvec4(7,1,3,5);
	glm::uvec4* polyNormalIndices = new glm::uvec4[6];
	polyNormalIndices[0] = glm::uvec4(1,2,3,4);
	polyNormalIndices[1] = glm::uvec4(5,6,7,8);
	polyNormalIndices[2] = glm::uvec4(9,10,11,12);
	polyNormalIndices[3] = glm::uvec4(13,14,15,16);
	polyNormalIndices[4] = glm::uvec4(17,18,19,20);
	polyNormalIndices[5] = glm::uvec4(21,22,23,24);
	glm::uvec4* polyUVIndices = new glm::uvec4[6];
	polyUVIndices[0] = glm::uvec4(1,2,3,4);
	polyUVIndices[1] = glm::uvec4(1,2,3,4);
	polyUVIndices[2] = glm::uvec4(3,4,1,2);
	polyUVIndices[3] = glm::uvec4(1,2,3,4);
	polyUVIndices[4] = glm::uvec4(1,2,3,4);
	polyUVIndices[5] = glm::uvec4(1,2,3,4);
	objCore::Obj* mesh = new objCore::Obj();
    mesh->m_numberOfVertices = 8;
    mesh->m_vertices = vertices;
    mesh->m_numberOfNormals = 24;
    mesh->m_normals = normals;
    mesh->m_numberOfUVs = 4;
    mesh->m_uvs = uvs;
    mesh->m_numberOfPolys = 6;
    mesh->m_polyVertexIndices = polyVertexIndices;
    mesh->m_polyNormalIndices = polyNormalIndices;
    mesh->m_polyUVIndices = polyUVIndices;
	return mesh;
}

objCore::Obj* Cube::Tesselate(const glm::vec3& lowerCorner, const glm::vec3& upperCorner){
	glm::vec3 scale = upperCorner-lowerCorner;
	glm::vec3 center = (upperCorner+lowerCorner)/2.0f;
	glm::mat4 transform = utilityCore::buildTransformationMatrix(center, glm::vec3(0,0,0), scale);

	objCore::Obj* o = Tesselate();
	unsigned int numberOfPoints = o->m_numberOfVertices;

	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,numberOfPoints),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				o->m_vertices[i] = glm::vec3(transform*glm::vec4(o->m_vertices[i], 1.0f));
			}
		}
	);

	return o;
}
}

