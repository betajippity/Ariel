// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: cube.cpp
// Implements cube.hpp

#include <tbb/tbb.h>
#include "cube.hpp"
#include "../utilities/utilities.h" 

using namespace geomCore;

//Default empty constructor defaults to 1 subdivs in axis and height
cube::cube(){	
}

//Boring empty destructor is boring
cube::~cube(){
}

//Literally just embed an obj because that's how we roll. Stupid but it works.
objCore::objContainer* cube::tesselate(){
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
	glm::vec4* polyVertexIndices = new glm::vec4[6];
	polyVertexIndices[0] = glm::vec4(1,2,4,3);
	polyVertexIndices[1] = glm::vec4(3,4,6,5);
	polyVertexIndices[2] = glm::vec4(5,6,8,7);
	polyVertexIndices[3] = glm::vec4(7,8,2,1);
	polyVertexIndices[4] = glm::vec4(2,8,6,4);
	polyVertexIndices[5] = glm::vec4(7,1,3,5);
	glm::vec4* polyNormalIndices = new glm::vec4[6];
	polyNormalIndices[0] = glm::vec4(1,2,3,4);
	polyNormalIndices[1] = glm::vec4(5,6,7,8);
	polyNormalIndices[2] = glm::vec4(9,10,11,12);
	polyNormalIndices[3] = glm::vec4(13,14,15,16);
	polyNormalIndices[4] = glm::vec4(17,18,19,20);
	polyNormalIndices[5] = glm::vec4(21,22,23,24);
	glm::vec4* polyUVIndices = new glm::vec4[6];
	polyUVIndices[0] = glm::vec4(1,2,3,4);
	polyUVIndices[1] = glm::vec4(1,2,3,4);
	polyUVIndices[2] = glm::vec4(3,4,1,2);
	polyUVIndices[3] = glm::vec4(1,2,3,4);
	polyUVIndices[4] = glm::vec4(1,2,3,4);
	polyUVIndices[5] = glm::vec4(1,2,3,4);
	objCore::obj* mesh = objCore::createObj(8, vertices, 24, normals, 4, uvs, 6, 
											polyVertexIndices, polyNormalIndices, polyUVIndices);
	objCore::objContainer* o = new objCore::objContainer(mesh);
	return o;
}

objCore::objContainer* cube::tesselate(const glm::vec3& lowerCorner, const glm::vec3& upperCorner){
	glm::vec3 scale = upperCorner-lowerCorner;
	glm::vec3 center = (upperCorner+lowerCorner)/2.0f;
	glm::mat4 transform = utilityCore::buildTransformationMatrix(center, glm::vec3(0,0,0), scale);

	objCore::objContainer* o = tesselate();
	unsigned int numberOfPoints = o->getObj()->numberOfVertices;

	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,numberOfPoints),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				o->getObj()->vertices[i] = glm::vec3(transform*glm::vec4(o->getObj()->vertices[i], 
													 1.0f));
			}
		}
	);

	return o;
}
