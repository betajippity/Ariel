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
#include <partio/Partio.h>

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

void scene::setPaths(const string& imagePath, const string& meshPath, const string& vdbPath,
					 const string& partioPath){
	this->imagePath = imagePath;
	this->meshPath = meshPath;
	this->vdbPath = vdbPath;
	this->partioPath = partioPath;
}

void scene::exportParticles(vector<fluidCore::particle*> particles, const float& maxd, 
							const int& frame, const bool& VDB, const bool& OBJ, 
							const bool& PARTIO){
	int particlesCount = particles.size();

	vector<fluidCore::particle*> sdfparticles;
	for(int i=0; i<particlesCount; i++){
		if(particles[i]->type==FLUID && !particles[i]->invalid){
			sdfparticles.push_back(particles[i]);
		}
	}
	int sdfparticlesCount = sdfparticles.size();
	
	string frameString = utilityCore::padString(4, utilityCore::convertIntToString(frame));

	if(PARTIO){
		string partiofilename = partioPath;
		vector<string> tokens = utilityCore::tokenizeString(partiofilename, ".");
		string ext = "." + tokens[tokens.size()-1];
		if(strcmp(ext.c_str(), ".gz")==0){
			ext = "." + tokens[tokens.size()-2] + ext;
		}

	    utilityCore::replaceString(partiofilename, ext, "."+frameString+ext);

		Partio::ParticlesDataMutable* partioData = Partio::create();
		partioData->addParticles(sdfparticlesCount);
		Partio::ParticleAttribute positionAttr = partioData->addAttribute("position", 
																		  Partio::VECTOR, 3);
		Partio::ParticleAttribute velocityAttr = partioData->addAttribute("v", Partio::VECTOR, 3);
		Partio::ParticleAttribute idAttr = partioData->addAttribute("id", Partio::INT, 1);

		for(int i=0; i<sdfparticlesCount; i++){
			float* pos = partioData->dataWrite<float>(positionAttr, i);
			pos[0] = sdfparticles[i]->p.x * maxd;
			pos[1] = sdfparticles[i]->p.y * maxd;
			pos[2] = sdfparticles[i]->p.z * maxd;
			float* vel = partioData->dataWrite<float>(velocityAttr, i);
			vel[0] = sdfparticles[i]->u.x;
			vel[1] = sdfparticles[i]->u.y;
			vel[2] = sdfparticles[i]->u.z;
			int* id = partioData->dataWrite<int>(idAttr, i);
			id[0] = i;			
		}

		Partio::write(partiofilename.c_str() ,*partioData);
		partioData->release();
	}

	if(VDB || OBJ){
		string vdbfilename = vdbPath;
	    utilityCore::replaceString(vdbfilename, ".vdb", "."+frameString+".vdb");

	    string objfilename = meshPath;
	    utilityCore::replaceString(objfilename, ".obj", "."+frameString+".obj");

		fluidCore::levelset* fluidSDF = new fluidCore::levelset(sdfparticles, maxd);

		if(VDB){
			fluidSDF->writeVDBGridToFile(vdbfilename);
		}

		if(OBJ){
			fluidSDF->writeObjToFile(objfilename);
		}
		delete fluidSDF;
	}

	
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

void scene::projectPointsToSolidSurface(vector<vec3>& points){
	vector<vec3> p1(points);
	solidLevelSet->projectPointsToSurface(p1);
	vector<vec3> p2(points);
	permaSolidLevelSet->projectPointsToSurface(p2);

	int pointsCount = points.size();

	for(int i=0; i<pointsCount; i++){
		float l1 = length(p1[i] - points[i]);
		float l2 = length(p2[i] - points[i]);
		if(l1<l2){
			points[i] = p1[i];
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

	// permaSolidLevelSet->writeVDBGridToFile("test.vdb");
}

// void scene::rebuildLiquidLevelSet(vector<fluidCore::particle*>& particles){
// 	delete liquidLevelSet;
// 	liquidLevelSet = new fluidCore::levelset(particles);
// }

void scene::generateParticles(vector<fluidCore::particle*>& particles, const vec3& dimensions, 
					   		  const float& density, fluidCore::particlegrid* pgrid, 
					   		  const int& frame){

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
							addParticle(vec3(x,y,z), FLUID, 3.0f/maxdimension, maxdimension, 
										particles, frame);
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
	                addParticle(vec3(x,y,z), SOLID, 3.0f/maxdimension, maxdimension, particles,
	                		    frame);
	            }
	        }
	    }
	}
    // cout << "Solid+Fluid particles: " << particles.size() << endl;
}

void scene::addParticle(const vec3& pos, const geomtype& type, const float& thickness, 
						const float& scale, vector<fluidCore::particle*>& particles, 
						const int& frame){
	bool inside = false;
	bool temp = false; //used to flag frame-variante solid particles

	if(type==FLUID){
		vec3 worldpos = pos*scale;
		if(liquidLevelSet->getInterpolatedCell(worldpos)<0.0f /*thickness*/){ 
			//TODO: figure out if we need this
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
			temp = true;
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
		p->temp = temp;
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
