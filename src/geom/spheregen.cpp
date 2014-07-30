// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: sphere.cpp
// Implements sphere.hpp

#include <tbb/tbb.h>
#include "spheregen.hpp"
#include "../utilities/utilities.h"

namespace geomCore{

//Default empty constructor defaults to 20 subdivs in axis and height
SphereGen::SphereGen(){	
	m_subdivs = 20;
}

//Constructor with options for presets
SphereGen::SphereGen(const unsigned int& subdivCount){
	m_subdivs = subdivCount;
}

//Boring empty destructor is boring
SphereGen::~SphereGen(){
}

void SphereGen::Tesselate(objCore::Obj* o, const glm::vec3& center, const float& radius){
	glm::vec3 scale = glm::vec3(radius*2.0f);
	glm::mat4 transform = utilityCore::buildTransformationMatrix(center, glm::vec3(0,0,0), scale);

	Tesselate(o);
	unsigned int numberOfPoints = o->m_numberOfVertices;

	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,numberOfPoints),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				o->m_vertices[i] = glm::vec3(transform*glm::vec4(o->m_vertices[i], 1.0f));
			}
		}
	);
}

/*Returns sphere mesh with specified subdiv counts. 
Axis and height must have a minimum of 3 subdivs and 
tesselate() will default to 3 if subdiv count is below 3.*/
//Yes, this function is total spaghetti and hacked together in many places. Will fix later. Maybe.
void SphereGen::Tesselate(objCore::Obj* o){
	unsigned int a = std::max(m_subdivs, 3);
	unsigned int h = std::max(m_subdivs, 3);
	unsigned int vertCount = (a*(h-1))+2;
	unsigned int faceCount = a*h;
	glm::vec3* vertices = new glm::vec3[vertCount];
	glm::vec3* normals = new glm::vec3[vertCount];
	glm::vec2* uvs = new glm::vec2[vertCount+(3*h)];
	glm::uvec4* polyVertexIndices = new glm::uvec4[faceCount];
	glm::uvec4* polyNormalIndices = new glm::uvec4[faceCount];
	glm::uvec4* polyUVIndices = new glm::uvec4[faceCount];
	//generate vertices and write into vertex array
	for(unsigned int x=1; x<h; x++){
		for(unsigned int y=0; y<a; y++){
			float angle1 = float(x)*float(1.0/a);
			float angle2 = float(y)*float(1.0/h);
			int i = ((x-1)*a) + y;
			vertices[i] = GetPointOnSphereByAngles(angle1,angle2);
		}
	}
	vertices[vertCount-2] = glm::vec3(0,-.5,0);
	vertices[vertCount-1] = glm::vec3(0,.5,0);
	//generate normals and uvs on sides of sphere
	float uvoffset = 0;
	for(unsigned int i=0; i<vertCount; i++){
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
	for(unsigned int i=0; i<vertCount; i++){
		uvs[i].x = uvs[i].x-uvoffset;
	}
	//generate wraparound uvs
	for(unsigned int i=0; i<h; i++){
		glm::vec2 uv;
		uv.x = 1;
		uv.y = uvs[i*h-1].y;
		uvs[vertCount+i] = uv;
	}
	//generate faces for sphere sides
	for(unsigned int x=1; x<h-1; x++){
		for(unsigned int y=0; y<a; y++){
			unsigned int i1 = ((x-1)*a) + y;
			unsigned int i2 = ((x-1)*a) + (y+1);
			unsigned int i3 = ((x)*a) + (y+1);
			unsigned int i4 = ((x)*a) + y;
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
	for(unsigned int x=0; x<h; x++){
		glm::uvec4 indices = polyVertexIndices[x];
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
		unsigned int uvindex = vertCount+h+x;
		glm::vec2 uv;
		uv.y=0;
		uv.x=uvs[(unsigned int)indices[2]].x+(.5/h);
		if(x==int(h/2)){
			uv.x = uv.x-(1.0f/h);
		}
		uvs[uvindex] = uv;
	}
	//generate faces and uvs for bottom pole
	for(unsigned int x=0; x<h; x++){
		unsigned int index = (h-1)*(a-2)-2+x;
		glm::uvec4 indices = polyVertexIndices[index];
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
		unsigned int uvindex = vertCount+(2*h)+x;
		glm::vec2 uv;
		uv.y=1;
		uv.x=uvs[(unsigned int)indices[0]].x+(.5/h);
		if(x==int(h/2)){
			uv.x = uv.x-(1.0f/h);
		}
		uvs[vertCount+(2*h)+x] = uv;
	}
	//flip uvs
	for(unsigned int i=0; i<vertCount+(3*h); i++){
		uvs[i].x = 1.0f-uvs[i].x;
	}
    o->m_numberOfVertices = vertCount;
    o->m_vertices = vertices;
    o->m_numberOfNormals = vertCount;
    o->m_normals = normals;
    o->m_numberOfUVs = vertCount+(3*h);
    o->m_uvs = uvs;
    o->m_numberOfPolys = faceCount;
    o->m_polyVertexIndices = polyVertexIndices;
    o->m_polyNormalIndices = polyNormalIndices;
    o->m_polyUVIndices = polyUVIndices;
}

glm::vec3 SphereGen::GetPointOnSphereByAngles(const float& angle1, const float& angle2){
	float x = sin(PI*angle1)*cos(TWO_PI*angle2);
	float y = sin(PI*angle1)*sin(TWO_PI*angle2);
	float z = cos(PI*angle1);
	return glm::normalize(glm::vec3(x,z,y))/2.0f;
}
}
