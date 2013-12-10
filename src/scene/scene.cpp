// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: scene.cpp
// Implements scene.hpp

#include "scene.hpp"

using namespace sceneCore;

scene::scene(){

}

scene::~scene(){

}

void scene::addGeom(geomCore::geom* object){
	if(object->getType()==SOLID){
		barriers.push_back(object);
	}else if(object->getType()==FLUID){
		barriers.push_back(object);
	}
}

void scene::generateParticles(vector<fluidCore::particle*>& particles, const vec3& dimensions, 
					   		  const float& density){
	
}