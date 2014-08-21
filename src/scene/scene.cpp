// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: scene.cpp
// Implements scene.hpp

#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>
#include <partio/Partio.h>
#include "scene.hpp"

namespace sceneCore{

Scene::Scene(){
	m_solidLevelSet = new fluidCore::LevelSet();
	m_liquidLevelSet = new fluidCore::LevelSet();
	m_highresSolidParticles = true;
	m_liquidParticleCount = 0;
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

void Scene::ExportParticles(std::vector<fluidCore::Particle*> particles, 
							const float& maxd, const int& frame, const bool& VDB, const bool& OBJ, 
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

std::vector<geomCore::Geom*>& Scene::GetSolidGeoms(){
	return m_solids;
}

std::vector<geomCore::Geom*>& Scene::GetLiquidGeoms(){
	return m_liquids;
}


void Scene::BuildLevelSets(const int& frame){
	delete m_liquidLevelSet;
	m_liquidLevelSet = new fluidCore::LevelSet();
	delete m_solidLevelSet;
	m_solidLevelSet = new fluidCore::LevelSet();

	unsigned int liquidObjectsCount = m_liquids.size();
	bool liquidSDFCreated = false;
	for(unsigned int i=0; i<liquidObjectsCount; i++){
		glm::mat4 transform;
		glm::mat4 inversetransform;
		if(m_liquids[i]->m_geom->GetTransforms((float)frame, transform, inversetransform)==true){
        	GeomType type = m_liquids[i]->m_geom->GetType();
        	if(type==MESH){
        		geomCore::MeshContainer* m = dynamic_cast<geomCore::MeshContainer*>
        									 			 (m_liquids[i]->m_geom);
        		objCore::Obj* o = &m->GetMeshFrame((float)frame)->m_basegeom;
        		if(liquidSDFCreated==false){
        			delete m_liquidLevelSet;
        			m_liquidLevelSet = new fluidCore::LevelSet(o, transform);
        			liquidSDFCreated = true;
        		}else{
        			fluidCore::LevelSet* objectSDF = new fluidCore::LevelSet(o, transform);
					m_liquidLevelSet->Merge(*objectSDF);
					delete objectSDF;
        		}
        	}else if(type==ANIMMESH){
        		geomCore::AnimatedMeshContainer* m=dynamic_cast<geomCore::AnimatedMeshContainer*>
        									 			 		  (m_liquids[i]->m_geom);
        		objCore::InterpolatedObj* o = &m->GetMeshFrame((float)frame)->m_basegeom;
        		float interpolationWeight = m->GetInterpolationWeight((float)frame);
        		if(liquidSDFCreated==false){
        			delete m_liquidLevelSet;
        			m_liquidLevelSet = new fluidCore::LevelSet(o, interpolationWeight,
        													   transform);
        			liquidSDFCreated = true;
        		}else{
        			fluidCore::LevelSet* objectSDF = new fluidCore::LevelSet(o,interpolationWeight,
        													   				 transform);
					m_liquidLevelSet->Merge(*objectSDF);
					delete objectSDF;
        		}
        	}
    	}
	}

	unsigned int solidObjectsCount = m_solids.size();
	bool solidSDFCreated = false;
	for(unsigned int i=0; i<solidObjectsCount; i++){
		glm::mat4 transform;
		glm::mat4 inversetransform;
		if(m_solids[i]->m_geom->GetTransforms((float)frame, transform, inversetransform)==true){
        	GeomType type = m_solids[i]->m_geom->GetType();
        	if(type==MESH){
        		geomCore::MeshContainer* m = dynamic_cast<geomCore::MeshContainer*>
        									 			 (m_solids[i]->m_geom);
        		objCore::Obj* o = &m->GetMeshFrame((float)frame)->m_basegeom;
        		if(solidSDFCreated==false){
        			delete m_solidLevelSet;
        			m_solidLevelSet = new fluidCore::LevelSet(o, transform);
        			solidSDFCreated = true;
        		}else{
        			fluidCore::LevelSet* objectSDF = new fluidCore::LevelSet(o, transform);
					m_solidLevelSet->Merge(*objectSDF);
					delete objectSDF;
        		}
        	}else if(type==ANIMMESH){
        		geomCore::AnimatedMeshContainer* m=dynamic_cast<geomCore::AnimatedMeshContainer*>
        									 			 		  (m_solids[i]->m_geom);
        		objCore::InterpolatedObj* o = &m->GetMeshFrame((float)frame)->m_basegeom;
        		float interpolationWeight = m->GetInterpolationWeight((float)frame);
        		if(solidSDFCreated==false){
        			delete m_solidLevelSet;
        			m_solidLevelSet = new fluidCore::LevelSet(o, interpolationWeight,
        													   transform);
        			solidSDFCreated = true;
        		}else{
        			fluidCore::LevelSet* objectSDF = new fluidCore::LevelSet(o,interpolationWeight,
        													   				 transform);
					m_solidLevelSet->Merge(*objectSDF);
					delete objectSDF;
        		}
        	}
    	}
	}
}

unsigned int Scene::GetLiquidParticleCount(){
	return m_liquidParticleCount;
}

void Scene::GenerateParticles(std::vector<fluidCore::Particle*>& particles,
							  const glm::vec3& dimensions, const float& density, 
							  fluidCore::ParticleGrid* pgrid, const int& frame){

	float maxdimension = glm::max(glm::max(dimensions.x, dimensions.y), dimensions.z);

	float thickness = 1.0f/maxdimension;
	float w = density*thickness;

    //store list of pointers to particles we need to delete for later deletion in the locked block
    std::vector<fluidCore::Particle*> particlesToDelete;
    particlesToDelete.reserve(m_solidParticles.size()+m_permaSolidParticles.size());
    particlesToDelete.insert(particlesToDelete.end(), m_solidParticles.begin(), m_solidParticles.end());
    particlesToDelete.insert(particlesToDelete.end(), m_permaSolidParticles.begin(),
                             m_permaSolidParticles.end());

    //swap-clear vectors
    tbb::concurrent_vector<fluidCore::Particle*>().swap(m_solidParticles);
    tbb::concurrent_vector<fluidCore::Particle*>().swap(m_permaSolidParticles);
    
    //place fluid particles
	if(m_liquids.size()>0){
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0,(dimensions.x+1)/density),
			[=](const tbb::blocked_range<unsigned int>& r){
				for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				  	for(unsigned int j = 0; j<(dimensions.y+1) / density; ++j){
						for(unsigned int k = 0; k<(dimensions.z+1) / density; ++k){
							float x = (i*w)+(w/2.0f);
							float y = (j*w)+(w/2.0f);
							float z = (k*w)+(w/2.0f);
							AddParticle(glm::vec3(x,y,z), FLUID, 3.0f/maxdimension, maxdimension, 
										frame);
						}
					}
				}
			}
		);
	}
	//std::cout << "Fluid particles: " << m_liquidParticles.size() << std::endl;

	if(m_highresSolidParticles==false){
		if(m_solids.size()>0){
		    w = 1.0f/maxdimension;
                tbb::parallel_for(tbb::blocked_range<unsigned int>(0,dimensions.x),
                [=](const tbb::blocked_range<unsigned int>& r){
                    for(unsigned int i=r.begin(); i!=r.end(); ++i){	
                        for(unsigned int j = 0; j < dimensions.y; j++) {
                            for(unsigned int k=0; k < dimensions.z; k++ ) {
                                float x = i*w+w/2.0f;
                                float y = j*w+w/2.0f;
                                float z = k*w+w/2.0f;
                                AddParticle(glm::vec3(x,y,z), SOLID, 3.0f/maxdimension,
                                            maxdimension, frame);
                            }
                        }
                    }
                }
			);
		}
	}else{
		if(m_solids.size()>0){
			tbb::parallel_for(tbb::blocked_range<unsigned int>(0,(dimensions.x+1)/density),
				[=](const tbb::blocked_range<unsigned int>& r){
					for(unsigned int i=r.begin(); i!=r.end(); ++i){	
					  	for(unsigned int j = 0; j<(dimensions.y+1) / density; ++j){
							for(unsigned int k = 0; k<(dimensions.z+1) / density; ++k){
								float x = (i*w)+(w/2.0f);
								float y = (j*w)+(w/2.0f);
								float z = (k*w)+(w/2.0f);
								AddParticle(glm::vec3(x,y,z), SOLID, 3.0f/maxdimension, 
											maxdimension, frame);
							}
						}
					}
				}
			);
		}
	}

	m_particleLock.lock();
   
    //delete old particles 
    unsigned int delParticleCount = particlesToDelete.size();
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,delParticleCount),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){
                delete particlesToDelete[i];
            }
        }
    );

    //add new particles to main particles list
	std::vector<fluidCore::Particle*>().swap(particles);
    particles.reserve(m_liquidParticles.size()+m_permaSolidParticles.size()+
                      m_solidParticles.size()); 
    particles.insert(particles.end(), m_liquidParticles.begin(), m_liquidParticles.end()); 
    particles.insert(particles.end(), m_permaSolidParticles.begin(), m_permaSolidParticles.end());
    particles.insert(particles.end(), m_solidParticles.begin(), m_solidParticles.end());
    m_liquidParticleCount = m_liquidParticles.size();

    //std::cout << "Solid+Fluid particles: " << particles.size() << std::endl;

    m_particleLock.unlock();
}

bool Scene::CheckPointInsideSolidGeom(const glm::vec3& p, const float& frame, 
									  unsigned int& solidGeomID){
	rayCore::Ray r;
	r.m_origin = p;
	r.m_frame = frame;
	r.m_direction = glm::normalize(glm::vec3(0,0,1));
	unsigned int solidGeomCount = m_solids.size();
	for(unsigned int i=0; i<solidGeomCount; i++){
		unsigned int hits = 0;
		spaceCore::HitCountTraverseAccumulator traverser(p);
		m_solids[i]->Intersect(r, traverser);
		bool hit = false;
		if(traverser.m_intersection.m_hit==true){
			if((traverser.m_numberOfHits)%2==1){
				solidGeomID = i;
				return true;	
			}
		}
	}
	return false;
}

bool Scene::CheckPointInsideLiquidGeom(const glm::vec3& p, const float& frame, 
									   unsigned int& liquidGeomID){
	rayCore::Ray r;
	r.m_origin = p;
	r.m_frame = frame;
	r.m_direction = glm::normalize(glm::vec3(0,0,1));
	unsigned int liquidGeomCount = m_liquids.size();
	for(unsigned int i=0; i<liquidGeomCount; i++){
		spaceCore::HitCountTraverseAccumulator traverser(p);
		m_liquids[i]->Intersect(r, traverser);
		bool hit = false;
		if(traverser.m_intersection.m_hit==true){
			if((traverser.m_numberOfHits)%2==1){
				liquidGeomID = i;
				return true;	
			}
		}
	}
	return false;
}

rayCore::Intersection Scene::IntersectSolidGeoms(const rayCore::Ray& r){
	rayCore::Intersection bestHit;
	unsigned int solidGeomCount = m_solids.size();
	for(unsigned int i=0; i<solidGeomCount; i++){
		spaceCore::TraverseAccumulator traverser;
		m_solids[i]->Intersect(r, traverser);
		bestHit = bestHit.CompareClosestAgainst(traverser.m_intersection, r.m_origin);
	}
	return bestHit;
}

void Scene::AddParticle(const glm::vec3& pos, const geomtype& type, const float& thickness, 
						const float& scale, const int& frame){
	bool inside = false;
	bool temp = false; //used to flag frame-variant solid particles

	if(type==FLUID){
		glm::vec3 worldpos = pos*scale;
		unsigned int liquidGeomID;
		if(CheckPointInsideLiquidGeom(pos*scale, frame, liquidGeomID)){
			inside = true;
			//if particles are in a solid, don't generate them
			unsigned int solidGeomID;
			if(CheckPointInsideSolidGeom(pos*scale, frame, solidGeomID)){
				inside = false;
			}
		}
	}else if(type==SOLID){
		unsigned int solidGeomID;
		if(CheckPointInsideSolidGeom(pos*scale, frame, solidGeomID)){
			inside = true;
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
		if(type==FLUID){
			m_liquidParticles.push_back(p);
		}else if(type==SOLID){
			if(p->m_temp==true){
				m_solidParticles.push_back(p);
			}else{
				m_permaSolidParticles.push_back(p);
			}
		}
	}
}

fluidCore::LevelSet* Scene::GetSolidLevelSet(){
	return m_solidLevelSet;	
}

fluidCore::LevelSet* Scene::GetLiquidLevelSet(){
	return m_liquidLevelSet;
}
}

