// Ariel: FLIP Fluid Simulator
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
	permaSolidLevelSet = new fluidCore::levelset();
	permaLiquidLevelSet = new fluidCore::levelset();
	permaLiquidSDFActive = false;
	permaSolidSDFActive = false;
}

scene::~scene(){
	delete solidLevelSet;
	delete liquidLevelSet;
}

void scene::setPaths(const string& imagePath, const string& meshPath, const string& vdbPath){
	this->imagePath = imagePath;
	this->meshPath = meshPath;
	this->vdbPath = vdbPath;
}

void scene::addSolidObject(objCore::objContainer* object, int startFrame, int endFrame){
	solidObjects.push_back(object);
	solidObjectFrameRanges.push_back(vec2(startFrame, endFrame));

	if(startFrame<0 && endFrame<0){
		if(permaSolidSDFActive==false){
			delete permaSolidLevelSet;
			permaSolidLevelSet = new fluidCore::levelset(object);
			permaSolidSDFActive = true;
		}else{
			fluidCore::levelset* objectSDF = new fluidCore::levelset(object);
			permaSolidLevelSet->merge(*objectSDF);
			delete objectSDF;
		}
	}
}

void scene::addLiquidObject(objCore::objContainer* object, int startFrame, int endFrame){
	liquidObjects.push_back(object);
	liquidObjectFrameRanges.push_back(vec2(startFrame, endFrame));
	
	if(startFrame<0 && endFrame<0){
		if(permaLiquidSDFActive==false){
			delete permaLiquidLevelSet;
			permaLiquidLevelSet = new fluidCore::levelset(object);
			permaLiquidSDFActive = true;
		}else{
			fluidCore::levelset* objectSDF = new fluidCore::levelset(object);
			permaLiquidLevelSet->merge(*objectSDF);
			delete objectSDF;
		}
	}
}

vector<objCore::objContainer*>& scene::getSolidObjects(){
	return solidObjects;
}

vector<objCore::objContainer*>& scene::getLiquidObjects(){
	return liquidObjects;
}

void scene::buildLevelSets(const int& frame){
	//first rebuild frame dependent levelsets, then merge with permanent sets if needed

	delete liquidLevelSet;
	liquidLevelSet = new fluidCore::levelset();
	delete solidLevelSet;
	solidLevelSet = new fluidCore::levelset();

	int liquidObjectsCount = liquidObjects.size();
	bool liquidSDFCreated = false;
	for(int i=0; i<liquidObjectsCount; i++){
		if( (frame<=liquidObjectFrameRanges[i][1] && frame>=liquidObjectFrameRanges[i][0]) ){
			if(liquidSDFCreated==false){
				delete liquidLevelSet;
				liquidLevelSet = new fluidCore::levelset(liquidObjects[i]);
				liquidSDFCreated = true;
			}else{
				fluidCore::levelset* objectSDF = new fluidCore::levelset(liquidObjects[i]);
				liquidLevelSet->merge(*objectSDF);
				delete objectSDF;
			}
		}
	}
	int solidObjectsCount = solidObjects.size();
	bool solidSDFCreated = false;
	for(int i=0; i<solidObjectsCount; i++){
		if( (frame<=solidObjectFrameRanges[i][1] && frame>=solidObjectFrameRanges[i][0]) ){
			if(solidSDFCreated==false){
				delete solidLevelSet;
				solidLevelSet = new fluidCore::levelset(solidObjects[i]);
				solidSDFCreated = true;
			}else{
				fluidCore::levelset* objectSDF = new fluidCore::levelset(solidObjects[i]);
				solidLevelSet->merge(*objectSDF);
				delete objectSDF;
			}
		}
	}

	if(permaLiquidSDFActive){
		if(!liquidSDFCreated){
			delete liquidLevelSet;
			liquidLevelSet = new fluidCore::levelset();
			liquidLevelSet->copy(*permaLiquidLevelSet);
		}else{
			liquidLevelSet->merge(*permaLiquidLevelSet);
		}
	}
}

// void scene::rebuildLiquidLevelSet(vector<fluidCore::particle*>& particles){
// 	delete liquidLevelSet;
// 	liquidLevelSet = new fluidCore::levelset(particles);
// }

void scene::generateParticles(vector<fluidCore::particle*>& particles, const vec3& dimensions, 
					   		  const float& density, fluidCore::particlegrid* pgrid, const int& frame){

	float maxdimension = glm::max(glm::max(dimensions.x, dimensions.y), dimensions.z);

	float thickness = 1.0f/maxdimension;
	float w = density*thickness;

	//place fluid particles
	if(liquidObjects.size()>0){
		for(int i=0; i<dimensions.x/density; i++){
			for(int j=0; j<dimensions.y/density; j++){
				for(int k=0; k<dimensions.z/density; k++){
					float x = (i*w)+(w/2.0f);
					float y = (j*w)+(w/2.0f);
					float z = (k*w)+(w/2.0f);


					if( x > thickness && x < 1.0-thickness &&
						y > thickness && y < 1.0-thickness &&
						z > thickness && z < 1.0-thickness ) {
							addParticle(vec3(x,y,z), FLUID, 3.0f/maxdimension, maxdimension, particles, 
										frame);
					}
				}
			}
		}
	}
	// cout << "Fluid particles: " << particles.size() << endl;

	if(solidObjects.size()>0){
	    w = 1.0f/maxdimension;
	    for( int i=0; i < dimensions.x; i++ ) {
	        for( int j=0; j < dimensions.y; j++ ) {
	            for( int k=0; k < dimensions.z; k++ ) {
	                float x = i*w+w/2.0f;
	                float y = j*w+w/2.0f;
	                float z = k*w+w/2.0f;
	                addParticle(vec3(x,y,z), SOLID, 3.0f/maxdimension, maxdimension, particles, frame);
	            }
	        }
	    }
	}
    // cout << "Solid+Fluid particles: " << particles.size() << endl;
}

void scene::addParticle(const vec3& pos, const geomtype& type, const float& thickness, const float& scale,
						vector<fluidCore::particle*>& particles, const int& frame){
	bool inside = false;

	if(type==FLUID){
		vec3 worldpos = pos*scale;
		if(liquidLevelSet->getInterpolatedCell(worldpos)<0.0f /*thickness*/){ //TODO: figure out if we need this
			inside = true;
		}
		//if particles are in a wall, don't generate them
		if(solidLevelSet->getInterpolatedCell(worldpos)<0.0f /*thickness*/){ 
			inside = false; 
		}	
		if(frame==0 && permaSolidSDFActive){
			if(permaSolidLevelSet->getInterpolatedCell(worldpos)<0.0f /*thickness*/){
				inside = false;
			}
		}

	}else if(type==SOLID){
		vec3 worldpos = pos*scale;
		if(solidLevelSet->getInterpolatedCell(worldpos)<0.0f /*thickness*/){
			inside = true;
		}	
		if(frame==0 && permaSolidSDFActive){
			if(permaSolidLevelSet->getInterpolatedCell(worldpos)<0.0f /*thickness*/){
				inside = true;
			}
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
		p->invalid = false;
		particles.push_back(p);
	}
}

fluidCore::levelset* scene::getSolidLevelSet(){
	return solidLevelSet;	
}

fluidCore::levelset* scene::getLiquidLevelSet(){
	return liquidLevelSet;
}

vec2 scene::getSolidFrameRange(const int& index){
	return solidObjectFrameRanges[index];
}

vec2 scene::getLiquidFrameRange(const int& index){
	return liquidObjectFrameRanges[index];
}
