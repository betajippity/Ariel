// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: sphere.cpp
// Implements sphere.hpp

#include <tbb/tbb.h>
#include "sphere.hpp"
#include "../utilities/utilities.h"

using namespace geomCore;

//Default empty constructor defaults to 20 subdivs in axis and height
sphere::sphere(){	
	subdivs = 20;
}

//Constructor with options for presets
sphere::sphere(int subdivCount){
	subdivs = subdivCount;
}

//Boring empty destructor is boring
sphere::~sphere(){
}

objCore::objContainer* sphere::tesselate(const glm::vec3& center, const float& radius){
	glm::vec3 scale = glm::vec3(radius*2.0f);
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

/*Returns sphere mesh with specified subdiv counts. 
Axis and height must have a minimum of 3 subdivs and 
tesselate() will default to 3 if subdiv count is below 3.*/
//Yes, this function is total spaghetti and hacked together in many places. Will fix later. Maybe.
objCore::objContainer* sphere::tesselate(){
	int a = std::max(subdivs, 3);
	int h = std::max(subdivs, 3);
	int vertCount = (a*(h-1))+2;
	int faceCount = a*h;
	glm::vec3* vertices = new glm::vec3[vertCount];
	glm::vec3* normals = new glm::vec3[vertCount];
	glm::vec2* uvs = new glm::vec2[vertCount+(3*h)];
	glm::vec4* polyVertexIndices = new glm::vec4[faceCount];
	glm::vec4* polyNormalIndices = new glm::vec4[faceCount];
	glm::vec4* polyUVIndices = new glm::vec4[faceCount];
	//generate vertices and write into vertex array
	for(int x=1; x<h; x++){
		for(int y=0; y<a; y++){
			float angle1 = float(x)*float(1.0/a);
			float angle2 = float(y)*float(1.0/h);
			int i = ((x-1)*a) + y;
			vertices[i] = getPointOnSphereByAngles(angle1,angle2);
		}
	}
	vertices[vertCount-2] = glm::vec3(0,-.5,0);
	vertices[vertCount-1] = glm::vec3(0,.5,0);
	//generate normals and uvs on sides of sphere
	float uvoffset = 0;
	for(int i=0; i<vertCount; i++){
		normals[i] = glm::normalize(vertices[i]);
		glm::vec2 uv;
		uv.x = 0.5 - (atan2(normals[i].z, normals[i].x)/TWO_PI) ;
		uv.y = 0.5 - (2.0*(asin(normals[i].y)/TWO_PI));
		if(i==int(a/2)){
			uvoffset = uv.x;
		}
		uvs[i] = uv;
	}
	//offset uvs
	for(int i=0; i<vertCount; i++){
		uvs[i].x = uvs[i].x-uvoffset;
	}
	//generate wraparound uvs
	for(int i=0; i<h; i++){
		glm::vec2 uv;
		uv.x = 1;
		uv.y = uvs[i*h-1].y;
		uvs[vertCount+i] = uv;
	}
	//generate faces for sphere sides
	for(int x=1; x<h-1; x++){
		for(int y=0; y<a; y++){
			int i1 = ((x-1)*a) + y;
			int i2 = ((x-1)*a) + (y+1);
			int i3 = ((x)*a) + (y+1);
			int i4 = ((x)*a) + y;
			//attach vertices at wraparound point
			if(y==a-1){
				i2 = ((x-1)*a) + (0);
				i3 = ((x)*a) + (0);
			}
			glm::vec4 indices = glm::vec4(i1+1,i2+1,i3+1,i4+1);
			polyVertexIndices[i1] = indices;
			polyNormalIndices[i1] = indices;
			//fix uvs at uv wraparound point
			if(y==a/2){
				indices[0] = vertCount+x+1;
				indices[3] = vertCount+x+2;
			}
			polyUVIndices[i1] = indices;
		}
	}
	//generate faces and uvs for top pole
	for(int x=0; x<h; x++){
		glm::vec4 indices = polyVertexIndices[x];
		indices[3] = -1;
		indices[2] = indices[0];
		indices[0] = vertCount;
		polyVertexIndices[faceCount-(a*2)+x] = indices;
		polyNormalIndices[faceCount-(a*2)+x] = indices;
		indices = polyUVIndices[x];
		indices[3] = -1;
		indices[2] = indices[0];
		indices[0] = vertCount+h+x+1;
		polyUVIndices[faceCount-(a*2)+x] = indices;
		int uvindex = vertCount+h+x;
		glm::vec2 uv;
		uv.y=0;
		uv.x=uvs[(int)indices[2]].x+(.5/h);
		if(x==int(h/2)){
			uv.x = uv.x-(1.0f/h);
		}
		uvs[uvindex] = uv;
	}
	//generate faces and uvs for bottom pole
	for(int x=0; x<h; x++){
		int index = (h-1)*(a-2)-2+x;
		glm::vec4 indices = polyVertexIndices[index];
		indices[1] = indices[2];
		indices[0] = indices[3];
		indices[3] = -1;
		indices[2] = vertCount-1;
		polyVertexIndices[faceCount-a+x] = indices;
		polyNormalIndices[faceCount-a+x] = indices;
		indices = polyUVIndices[vertCount-(h*2)-2+x];
		glm::vec4 i2;
		i2[0] = indices[3];
		i2[1] = indices[2];
		i2[2] = vertCount+(h*2)+x+1;
		i2[3] = -1;
		polyUVIndices[faceCount-a+x] = i2;
		int uvindex = vertCount+(2*h)+x;
		glm::vec2 uv;
		uv.y=1;
		uv.x=uvs[(int)indices[0]].x+(.5/h);
		if(x==int(h/2)){
			uv.x = uv.x-(1.0f/h);
		}
		uvs[vertCount+(2*h)+x] = uv;
	}
	//flip uvs
	for(int i=0; i<vertCount+(3*h); i++){
		uvs[i].x = 1.0f-uvs[i].x;
	}
	objCore::obj* mesh = objCore::createObj(vertCount, vertices, vertCount, normals, 
											vertCount+(3*h), uvs, faceCount, polyVertexIndices,
											polyNormalIndices, polyUVIndices);
	objCore::objContainer* o = new objCore::objContainer(mesh);
	return o;
}

glm::vec3 sphere::getPointOnSphereByAngles(float angle1, float angle2){
	float x = sin(PI*angle1)*cos(TWO_PI*angle2);
	float y = sin(PI*angle1)*sin(TWO_PI*angle2);
	float z = cos(PI*angle1);
	return glm::normalize(glm::vec3(x,z,y))/2.0f;
}
