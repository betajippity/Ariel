// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: scene.cpp
// Implements scene.hpp

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Composite.h>
#include <partio/Partio.h>
#include "scene.hpp"

namespace sceneCore{

Scene::Scene(){
	m_solidLevelSet = new fluidCore::LevelSet();
	m_liquidLevelSet = new fluidCore::LevelSet();
	m_permaSolidLevelSet = new fluidCore::LevelSet();
	m_permaLiquidLevelSet = new fluidCore::LevelSet();
	m_permaLiquidSDFActive = false;
	m_permaSolidSDFActive = false;
}

Scene::~Scene(){
	delete m_solidLevelSet;
	delete m_liquidLevelSet;
}

void Scene::SetPaths(const std::string& imagePath, const std::string& meshPath, 
					 const std::string& vdbPath, const std::string& partioPath){
	m_imagePath = imagePath;
	m_meshPath = meshPath;
	m_vdbPath = vdbPath;
	m_partioPath = partioPath;
}

void Scene::ExportParticles(std::vector<fluidCore::Particle*> particles, const float& maxd,
							const int& frame, const bool& VDB, const bool& OBJ, 
							const bool& PARTIO){
	unsigned int particlesCount = particles.size();

	std::vector<fluidCore::Particle*> sdfparticles;
	for(unsigned int i = 0; i<particlesCount; i++){
		if(particles[i]->m_type==FLUID && !particles[i]->m_invalid){
			sdfparticles.push_back(particles[i]);
		}
	}
	int sdfparticlesCount = sdfparticles.size();
	
	std::string frameString = utilityCore::padString(4, utilityCore::convertIntToString(frame));

	if(PARTIO){
		std::string partiofilename = m_partioPath;
		std::vector<std::string> tokens = utilityCore::tokenizeString(partiofilename, ".");
		std::string ext = "." + tokens[tokens.size()-1];
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

		for(unsigned int i = 0; i<sdfparticlesCount; i++){
			float* pos = partioData->dataWrite<float>(positionAttr, i);
			pos[0] = sdfparticles[i]->m_p.x * maxd;
			pos[1] = sdfparticles[i]->m_p.y * maxd;
			pos[2] = sdfparticles[i]->m_p.z * maxd;
			float* vel = partioData->dataWrite<float>(velocityAttr, i);
			vel[0] = sdfparticles[i]->m_u.x;
			vel[1] = sdfparticles[i]->m_u.y;
			vel[2] = sdfparticles[i]->m_u.z;
			int* id = partioData->dataWrite<int>(idAttr, i);
			id[0] = i;			
		}

		Partio::write(partiofilename.c_str() ,*partioData);
		partioData->release();
	}

	if(VDB || OBJ){
		std::string vdbfilename = m_vdbPath;
	    utilityCore::replaceString(vdbfilename, ".vdb", "."+frameString+".vdb");

	    std::string objfilename = m_meshPath;
	    utilityCore::replaceString(objfilename, ".obj", "."+frameString+".obj");

		fluidCore::LevelSet* fluidSDF = new fluidCore::LevelSet(sdfparticles, maxd);

		if(VDB){
			fluidSDF->WriteVDBGridToFile(vdbfilename);
		}

		if(OBJ){
			fluidSDF->WriteObjToFile(objfilename);
		}
		delete fluidSDF;
	}

	
}

void Scene::AddExternalForce(glm::vec3 force){
	m_externalForces.push_back(force);
}

std::vector<glm::vec3>& Scene::GetExternalForces(){
	return m_externalForces;
}

void Scene::AddSolidObject(objCore::Obj* object, const int& startFrame, const int& endFrame){
	m_solidObjects.push_back(object);
	m_solidObjectFrameRanges.push_back(glm::vec2(startFrame, endFrame));

	if(startFrame<0 && endFrame<0){
		if(m_permaSolidSDFActive==false){
			delete m_permaSolidLevelSet;
			m_permaSolidLevelSet = new fluidCore::LevelSet(object);
			m_permaSolidSDFActive = true;
		}else{
			fluidCore::LevelSet* objectSDF = new fluidCore::LevelSet(object);
			m_permaSolidLevelSet->Merge(*objectSDF);
			delete objectSDF;
		}
	}
}

void Scene::AddLiquidObject(objCore::Obj* object, const int& startFrame, const int& endFrame){
	m_liquidObjects.push_back(object);
	m_liquidObjectFrameRanges.push_back(glm::vec2(startFrame, endFrame));
	
	if(startFrame<0 && endFrame<0){
		if(m_permaLiquidSDFActive==false){
			delete m_permaLiquidLevelSet;
			m_permaLiquidLevelSet = new fluidCore::LevelSet(object);
			m_permaLiquidSDFActive = true;
		}else{
			fluidCore::LevelSet* objectSDF = new fluidCore::LevelSet(object);
			m_permaLiquidLevelSet->Merge(*objectSDF);
			delete objectSDF;
		}
	}
}

void Scene::ProjectPointsToSolidSurface(std::vector<glm::vec3>& points){
	std::vector<glm::vec3> p1(points);
	m_solidLevelSet->ProjectPointsToSurface(p1);
	std::vector<glm::vec3> p2(points);
	m_permaSolidLevelSet->ProjectPointsToSurface(p2);

	unsigned int pointsCount = points.size();

	for(unsigned int i = 0; i<pointsCount; i++){
		float l1 = glm::length(p1[i] - points[i]);
		float l2 = glm::length(p2[i] - points[i]);
		if(l1<l2){
			points[i] = p1[i];
		}
	}
}

std::vector<objCore::Obj*>& Scene::GetSolidObjects(){
	return m_solidObjects;
}

std::vector<objCore::Obj*>& Scene::GetLiquidObjects(){
	return m_liquidObjects;
}

void Scene::BuildLevelSets(const int& frame){
	//first rebuild frame dependent levelsets, then merge with permanent sets if needed

	delete m_liquidLevelSet;
	m_liquidLevelSet = new fluidCore::LevelSet();
	delete m_solidLevelSet;
	m_solidLevelSet = new fluidCore::LevelSet();

	unsigned int liquidObjectsCount = m_liquidObjects.size();
	bool liquidSDFCreated = false;
	for(unsigned int i = 0; i<liquidObjectsCount; i++){
		if( (frame<=m_liquidObjectFrameRanges[i][1] && frame>=m_liquidObjectFrameRanges[i][0]) ){
			if(liquidSDFCreated==false){
				delete m_liquidLevelSet;
				m_liquidLevelSet = new fluidCore::LevelSet(m_liquidObjects[i]);
				liquidSDFCreated = true;
			}else{
				fluidCore::LevelSet* objectSDF = new fluidCore::LevelSet(m_liquidObjects[i]);
				m_liquidLevelSet->Merge(*objectSDF);
				delete objectSDF;
			}
		}
	}
	unsigned int solidObjectsCount = m_solidObjects.size();
	bool solidSDFCreated = false;
	for (unsigned int i = 0; i<solidObjectsCount; i++){
		if( (frame<=m_solidObjectFrameRanges[i][1] && frame>=m_solidObjectFrameRanges[i][0]) ){
			if(solidSDFCreated==false){
				delete m_solidLevelSet;
				m_solidLevelSet = new fluidCore::LevelSet(m_solidObjects[i]);
				solidSDFCreated = true;
			}else{
				fluidCore::LevelSet* objectSDF = new fluidCore::LevelSet(m_solidObjects[i]);
				m_solidLevelSet->Merge(*objectSDF);
				delete objectSDF;
			}
		}
	}

	if(m_permaLiquidSDFActive){
		if(!liquidSDFCreated){
			delete m_liquidLevelSet;
			m_liquidLevelSet = new fluidCore::LevelSet();
			m_liquidLevelSet->Copy(*m_permaLiquidLevelSet);
		}else{
			m_liquidLevelSet->Merge(*m_permaLiquidLevelSet);
		}
	}

	// permaSolidLevelSet->writeVDBGridToFile("test.vdb");
}

// void scene::rebuildLiquidLevelSet(std::vector<fluidCore::particle*>& particles){
// 	delete liquidLevelSet;
// 	liquidLevelSet = new fluidCore::levelset(particles);
// }

void Scene::GenerateParticles(std::vector<fluidCore::Particle*>& particles, 
							  const glm::vec3& dimensions, const float& density, 
							  fluidCore::ParticleGrid* pgrid, const int& frame){

	float maxdimension = glm::max(glm::max(dimensions.x, dimensions.y), dimensions.z);

	float thickness = 1.0f/maxdimension;
	float w = density*thickness;

	//place fluid particles
	if(m_liquidObjects.size()>0){
		for(unsigned int i = 0; i<dimensions.x / density; i++){
			for(unsigned int j = 0; j<dimensions.y / density; j++){
				for(unsigned int k = 0; k<dimensions.z / density; k++){
					float x = (i*w)+(w/2.0f);
					float y = (j*w)+(w/2.0f);
					float z = (k*w)+(w/2.0f);


					if( x > thickness && x < 1.0-thickness &&
						y > thickness && y < 1.0-thickness &&
						z > thickness && z < 1.0-thickness ) {
							AddParticle(glm::vec3(x,y,z), FLUID, 3.0f/maxdimension, maxdimension, 
										particles, frame);
					}
				}
			}
		}
	}
	// std::cout << "Fluid particles: " << particles.size() << std::endl;

	if(m_solidObjects.size()>0){
	    w = 1.0f/maxdimension;
		for(unsigned int i = 0; i < dimensions.x; i++) {
			for(unsigned int j = 0; j < dimensions.y; j++) {
	            for(unsigned int k=0; k < dimensions.z; k++ ) {
	                float x = i*w+w/2.0f;
	                float y = j*w+w/2.0f;
	                float z = k*w+w/2.0f;
	                AddParticle(glm::vec3(x,y,z), SOLID, 3.0f/maxdimension, maxdimension, 
	                			particles, frame);
	            }
	        }
	    }
	}
    // std::cout << "Solid+Fluid particles: " << particles.size() << std::endl;
}

void Scene::AddParticle(const glm::vec3& pos, const geomtype& type, const float& thickness, 
						const float& scale, std::vector<fluidCore::Particle*>& particles, 
						const int& frame){
	bool inside = false;
	bool temp = false; //used to flag frame-variante solid particles

	if(type==FLUID){
		glm::vec3 worldpos = pos*scale;
		if(m_liquidLevelSet->GetInterpolatedCell(worldpos)<0.0f /*thickness*/){ 
			//TODO: figure out if we need this
			inside = true;
		}
		//if particles are in a wall, don't generate them
		if(m_solidLevelSet->GetInterpolatedCell(worldpos)<0.0f /*thickness*/){ 
			inside = false; 
		}	
		if(frame==0 && m_permaSolidSDFActive){
			if(m_permaSolidLevelSet->GetInterpolatedCell(worldpos)<0.0f /*thickness*/){
				inside = false;
			}
		}

	}else if(type==SOLID){
		glm::vec3 worldpos = pos*scale;
		if(m_solidLevelSet->GetInterpolatedCell(worldpos)<0.0f /*thickness*/){
			inside = true;
			temp = true;
		}	
		if(frame==0 && m_permaSolidSDFActive){
			if(m_permaSolidLevelSet->GetInterpolatedCell(worldpos)<0.0f /*thickness*/){
				inside = true;
			}
		}
	}

	if(inside){
		fluidCore::Particle* p = new fluidCore::Particle;
		p->m_p = pos;
		p->m_u = glm::vec3(0,0,0);
		p->m_n = glm::vec3(0,0,0);
		p->m_density = 10.0f;
		p->m_type = type;
		p->m_mass = 1.0f;
		p->m_invalid = false;
		p->m_temp = temp;
		particles.push_back(p);
	}
}

fluidCore::LevelSet* Scene::GetSolidLevelSet(){
	return m_solidLevelSet;	
}

fluidCore::LevelSet* Scene::GetLiquidLevelSet(){
	return m_liquidLevelSet;
}

glm::vec2 Scene::GetSolidFrameRange(const int& index){
	return m_solidObjectFrameRanges[index];
}

glm::vec2 Scene::GetLiquidFrameRange(const int& index){
	return m_liquidObjectFrameRanges[index];
}
}
