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
	objects.push_back(object);
}

void scene::generateParticles(vector<fluidCore::particle*>& particles, const vec3& dimensions, 
					   		  const float& density){

	float maxdimension = glm::max(glm::max(dimensions.x, dimensions.y), dimensions.z);

	float thickness = 1.0f/maxdimension;

	float w = density*thickness;
	for(int i=0; i<dimensions.x/density; i++){
		for(int j=0; j<dimensions.y/density; j++){
			for(int k=0; k<dimensions.z/density; k++){
				float x = i*w+w/2.0f;
				float y = j*w+w/2.0f;
				float z = k*w+w/2.0f;

				if( x > thickness && x < 1.0-thickness &&
					y > thickness && y < 1.0-thickness &&
					z > thickness && z < 1.0-thickness ) {
						addParticle(vec3(x,y,z), FLUID, 3.0f/maxdimension, 1.0f/maxdimension, particles);
				}
			}
		}
	}
	cout << particles.size() << endl;
}

void scene::addParticle(const vec3& pos, const geomtype& type, const float& thickness, const float& scale,
						vector<fluidCore::particle*>& particles){
	geomCore::geom* boundingObject = NULL;
	for(int n=0; n<objects.size(); n++){

		bool found = false;
		
		if(objects[n]->isPointInside(pos, scale)==true){
			found = true;
			if(objects[n]->isPointInsideWithThickness(pos, thickness, scale)==true){
				boundingObject = NULL;
			}
		}

		if(found){
			if(objects[n]->getType()==type){
				boundingObject = objects[n];
				break;
			}
		}
	}

	if(boundingObject){
		fluidCore::particle* p = new fluidCore::particle;
		p->p = pos;
		p->u = vec3(0,0,0);
		p->n = vec3(0,0,0);
		p->density = 10.0f;
		p->type = type;
		p->mass = 1.0f;
		particles.push_back(p);
	}

}
