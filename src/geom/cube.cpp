// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: cube.cp p
// Impleme nts cube.hpp

#include "cube.hpp"
#include "../utilities/utilities.h" 
#include <omp.h>

using namespace geomCore;

//Default empty constructor defaults to 1 subdivs in axis and height
cube::cube(){	
}

//Boring empty destructor is boring
cube::~cube(){
}

//Literally just embed an obj because that's how we roll. Stupid but it works.
objCore::objContainer* cube::tesselate(){
	vec3* vertices = new vec3[8];
	vertices[0] = vec3(-0.5, -0.5,  0.5);
	vertices[1] = vec3( 0.5, -0.5,  0.5);
	vertices[2] = vec3(-0.5,  0.5,  0.5);
	vertices[3] = vec3( 0.5,  0.5,  0.5);
	vertices[4] = vec3(-0.5,  0.5, -0.5);
	vertices[5] = vec3( 0.5,  0.5, -0.5);
	vertices[6] = vec3(-0.5, -0.5, -0.5);
	vertices[7] = vec3( 0.5, -0.5, -0.5);
	vec3* normals = new vec3[24];
	normals[ 0] = vec3( 0.000000,  0.000000,  1.000000);
	normals[ 1] = vec3( 0.000000,  0.000000,  1.000000);
	normals[ 2] = vec3( 0.000000,  0.000000,  1.000000);
	normals[ 3] = vec3( 0.000000,  0.000000,  1.000000);
	normals[ 4] = vec3( 0.000000,  1.000000,  0.000000);
	normals[ 5] = vec3( 0.000000,  1.000000,  0.000000);
	normals[ 6] = vec3( 0.000000,  1.000000,  0.000000);
	normals[ 7] = vec3( 0.000000,  1.000000,  0.000000);
	normals[ 8] = vec3( 0.000000,  0.000000, -1.000000);
	normals[ 9] = vec3( 0.000000,  0.000000, -1.000000);
	normals[10] = vec3( 0.000000,  0.000000, -1.000000);
	normals[11] = vec3( 0.000000,  0.000000, -1.000000);
	normals[12] = vec3( 0.000000, -1.000000,  0.000000);
	normals[13] = vec3( 0.000000, -1.000000,  0.000000);
	normals[14] = vec3( 0.000000, -1.000000,  0.000000);
	normals[15] = vec3( 0.000000, -1.000000,  0.000000);
	normals[16] = vec3( 1.000000,  0.000000,  0.000000);
	normals[17] = vec3( 1.000000,  0.000000,  0.000000);
	normals[18] = vec3( 1.000000,  0.000000,  0.000000);
	normals[19] = vec3( 1.000000,  0.000000,  0.000000);
	normals[20] = vec3(-1.000000,  0.000000,  0.000000);
	normals[21] = vec3(-1.000000,  0.000000,  0.000000);
	normals[22] = vec3(-1.000000,  0.000000,  0.000000);
	normals[23] = vec3(-1.000000,  0.000000,  0.000000);
	vec2* uvs = new vec2[4];
	uvs[0] = vec2( 0, 0);
	uvs[1] = vec2( 1, 0);
	uvs[2] = vec2( 1, 1);
	uvs[3] = vec2( 0, 1);
	vec4* polyVertexIndices = new vec4[6];
	polyVertexIndices[0] = vec4(1,2,4,3);
	polyVertexIndices[1] = vec4(3,4,6,5);
	polyVertexIndices[2] = vec4(5,6,8,7);
	polyVertexIndices[3] = vec4(7,8,2,1);
	polyVertexIndices[4] = vec4(2,8,6,4);
	polyVertexIndices[5] = vec4(7,1,3,5);
	vec4* polyNormalIndices = new vec4[6];
	polyNormalIndices[0] = vec4(1,2,3,4);
	polyNormalIndices[1] = vec4(5,6,7,8);
	polyNormalIndices[2] = vec4(9,10,11,12);
	polyNormalIndices[3] = vec4(13,14,15,16);
	polyNormalIndices[4] = vec4(17,18,19,20);
	polyNormalIndices[5] = vec4(21,22,23,24);
	vec4* polyUVIndices = new vec4[6];
	polyUVIndices[0] = vec4(1,2,3,4);
	polyUVIndices[1] = vec4(1,2,3,4);
	polyUVIndices[2] = vec4(3,4,1,2);
	polyUVIndices[3] = vec4(1,2,3,4);
	polyUVIndices[4] = vec4(1,2,3,4);
	polyUVIndices[5] = vec4(1,2,3,4);
	objCore::obj* mesh = objCore::createObj(8, vertices, 24, normals, 4, uvs, 6, 
											polyVertexIndices, polyNormalIndices, polyUVIndices);
	objCore::objContainer* o = new objCore::objContainer(mesh);
	return o;
}

objCore::objContainer* cube::tesselate(const vec3& lowerCorner, const vec3& upperCorner){
	vec3 scale = upperCorner-lowerCorner;
	vec3 center = (upperCorner+lowerCorner)/2.0f;
	mat4 transform = utilityCore::buildTransformationMatrix(center, vec3(0,0,0), scale);

	objCore::objContainer* o = tesselate();
	int numberOfPoints = o->getObj()->numberOfVertices;

	#pragma omp parallel for
	for(int i=0; i<numberOfPoints; i++){
		o->getObj()->vertices[i] = vec3(transform*vec4(o->getObj()->vertices[i], 1.0f));
	}

	return o;
}
