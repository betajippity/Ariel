// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: scene.cpp
// Implements scene.hpp

#include "scene.hpp"
#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Composite.h>

using namespace sceneCore;

scene::scene(){
	solidLevelSet = new fluidCore::levelset();
	liquidLevelSet = new fluidCore::levelset();
}

scene::~scene(){
	delete solidLevelSet;
	delete liquidLevelSet;
}

void scene::addSolidObject(objCore::objContainer* object){
	solidObjects.push_back(object);
	if(solidObjects.size()==1){
		delete solidLevelSet;
		solidLevelSet = new fluidCore::levelset(object);
	}else{
		fluidCore::levelset* objectSDF = new fluidCore::levelset(object);
		solidLevelSet->merge(*objectSDF);
		delete objectSDF;
	}
}

void scene::addLiquidObject(objCore::objContainer* object){
	liquidObjects.push_back(object);
	if(liquidObjects.size()==1){
		delete liquidLevelSet;
		liquidLevelSet = new fluidCore::levelset(object);
	}else{
		fluidCore::levelset* objectSDF = new fluidCore::levelset(object);
		liquidLevelSet->merge(*objectSDF);
		delete objectSDF;
	}
}

vector<objCore::objContainer*>& scene::getSolidObjects(){
	return solidObjects;
}

vector<objCore::objContainer*>& scene::getLiquidObjects(){
	return liquidObjects;
}

void scene::rebuildLiquidLevelSet(vector<fluidCore::particle*>& particles){
	delete liquidLevelSet;
	liquidLevelSet = new fluidCore::levelset(particles);
}

void scene::generateParticles(vector<fluidCore::particle*>& particles, const vec3& dimensions, 
					   		  const float& density, fluidCore::particlegrid* pgrid){

	float maxdimension = glm::max(glm::max(dimensions.x, dimensions.y), dimensions.z);

	float thickness = 1.0f/maxdimension;

	//place fluid particles
	float w = density*thickness;
	for(int i=0; i<dimensions.x/density; i++){
		for(int j=0; j<dimensions.y/density; j++){
			for(int k=0; k<dimensions.z/density; k++){
				float x = (i*w)+(w/2.0f);
				float y = (j*w)+(w/2.0f);
				float z = (k*w)+(w/2.0f);
				if( x > thickness && x < 1.0-thickness &&
					y > thickness && y < 1.0-thickness &&
					z > thickness && z < 1.0-thickness ) {
						addParticle(vec3(x,y,z), FLUID, 3.0f/maxdimension, maxdimension, particles);
				}
			}
		}
	}
	cout << "Fluid particles: " << particles.size() << endl;

	// for( vector<particle *>::iterator iter=particles.begin(); iter!=particles.end(); ){
	// 	particle &p = **iter;

}

void scene::addParticle(const vec3& pos, const geomtype& type, const float& thickness, const float& scale,
						vector<fluidCore::particle*>& particles){
	bool inside = false;

	if(type==FLUID){
		vec3 worldpos = pos*scale;
		if(liquidLevelSet->getInterpolatedCell(worldpos)<thickness){
			inside = true;
		}
		//if particles are in a wall, don't generate them
		if(solidLevelSet->getInterpolatedCell(worldpos)<thickness){ 
			inside = false; 
		}	

	}else if(type==SOLID){
		vec3 worldpos = pos*scale;
		if(solidLevelSet->getInterpolatedCell(worldpos)<thickness){
			inside = true;
		}	
	}

	if(inside){
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

fluidCore::levelset* scene::getSolidLevelSet(){
	return solidLevelSet;	
}

fluidCore::levelset* scene::getLiquidLevelSet(){
	return liquidLevelSet;
}

