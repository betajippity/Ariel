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
    particlesToDelete.reserve(m_solidParticles.size());
    particlesToDelete.insert(particlesToDelete.end(), m_solidParticles.begin(), 
    						 m_solidParticles.end());

    tbb::concurrent_vector<fluidCore::Particle*>().swap(m_solidParticles);
    
    //place fluid particles
    //for each fluid geom in the frame, loop through voxels in the geom's AABB to place particles
	unsigned int liquidCount = m_liquids.size();
    for(unsigned int l=0; l<liquidCount; ++l){

        spaceCore::Aabb liquidaabb = m_liquids[l]->m_geom->GetAabb(frame);   
        glm::vec3 liquidvelocity = m_liquidStartingVelocities[l]; 

        if(m_liquids[l]->m_geom->IsInFrame(frame)){
            //clip AABB to sim boundaries, account for density
            glm::vec3 lmin = glm::floor(liquidaabb.m_min);
            glm::vec3 lmax = glm::ceil(liquidaabb.m_max);
            lmin = glm::max(lmin, glm::vec3(0.0f))/density;             
            lmax = glm::min(lmax, dimensions+glm::vec3(1.0f))/density;
		    //place particles in AABB
            tbb::parallel_for(tbb::blocked_range<unsigned int>(lmin.x,lmax.x),
			    [=](const tbb::blocked_range<unsigned int>& r){
				    for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				  	    for(unsigned int j = lmin.y; j<lmax.y; ++j){
						    for(unsigned int k = lmin.z; k<lmax.z; ++k){
							    float x = (i*w)+(w/2.0f);
							    float y = (j*w)+(w/2.0f);
							    float z = (k*w)+(w/2.0f);
							    AddLiquidParticle(glm::vec3(x,y,z), liquidvelocity, 
							    				  3.0f/maxdimension, maxdimension, 
							    				  frame, m_liquids[l]->m_id);
                            }
                        }
                    }
                }
            );
        }   
    }
	unsigned int solidCount = m_solids.size();
    for(unsigned int l=0; l<solidCount; ++l){
        spaceCore::Aabb solidaabb = m_solids[l]->m_geom->GetAabb(frame);        
        if((frame==0 && m_solids[l]->m_geom->IsDynamic()==false) || 
    	   (m_solids[l]->m_geom->IsDynamic()==true && m_solids[l]->m_geom->IsInFrame(frame))){
            //clip AABB to sim boundaries, account for density
            glm::vec3 lmin = glm::floor(solidaabb.m_min);
            glm::vec3 lmax = glm::ceil(solidaabb.m_max);
            lmin = glm::max(lmin, glm::vec3(0.0f))/density;             
            lmax = glm::min(lmax, dimensions+glm::vec3(1.0f))/density;
		    //place particles in AABB
            tbb::parallel_for(tbb::blocked_range<unsigned int>(lmin.x,lmax.x),
			    [=](const tbb::blocked_range<unsigned int>& r){
				    for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				  	    for(unsigned int j = lmin.y; j<lmax.y; ++j){
						    for(unsigned int k = lmin.z; k<lmax.z; ++k){
							    float x = (i*w)+(w/2.0f);
							    float y = (j*w)+(w/2.0f);
							    float z = (k*w)+(w/2.0f);
							    AddSolidParticle(glm::vec3(x,y,z), 3.0f/maxdimension, 
                                            	 maxdimension, frame, m_solids[l]->m_id);
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

bool Scene::CheckPointInsideGeomByID(const glm::vec3& p, const float& frame, 
									 const unsigned int& geomID){
	if(geomID<m_geoms.size()){
		rayCore::Ray r;
		r.m_origin = p;
		r.m_frame = frame;
		r.m_direction = glm::normalize(glm::vec3(0,0,1));
		unsigned int hits = 0;
		spaceCore::HitCountTraverseAccumulator traverser(p);
		m_geoms[geomID].Intersect(r, traverser);
		bool hit = false;
		if(traverser.m_intersection.m_hit==true){
			if((traverser.m_numberOfHits)%2==1){
				return true;	
			}
		}
	}
	return false;
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

void Scene::AddLiquidParticle(const glm::vec3& pos, const glm::vec3& vel, const float& thickness, 
							  const float& scale, const int& frame, 
							  const unsigned int& liquidGeomID){
	glm::vec3 worldpos = pos*scale;
	if(CheckPointInsideGeomByID(worldpos, frame, liquidGeomID)==true){
		//if particles are in a solid, don't generate them
		unsigned int solidGeomID;
		if(CheckPointInsideSolidGeom(worldpos, frame, solidGeomID)==false){
			fluidCore::Particle* p = new fluidCore::Particle;
			p->m_p = pos;
			p->m_u = vel;
			p->m_n = glm::vec3(0.0f);
			p->m_density = 10.0f;
			p->m_type = FLUID;
			p->m_mass = 1.0f;
			p->m_invalid = false;
			m_liquidParticles.push_back(p);
		}
	}
}

void Scene::AddSolidParticle(const glm::vec3& pos, const float& thickness, const float& scale, 
							 const int& frame, const unsigned int& solidGeomID){
	glm::vec3 worldpos = pos*scale;
	if(CheckPointInsideGeomByID(worldpos, frame, solidGeomID)==true){
		fluidCore::Particle* p = new fluidCore::Particle;
		p->m_p = pos;
		p->m_u = glm::vec3(0.0f);
		p->m_n = glm::vec3(0.0f);
		p->m_density = 10.0f;
		p->m_type = SOLID;
		p->m_mass = 10.0f;
		p->m_invalid = false;
		if(frame==0 && m_geoms[solidGeomID].m_geom->IsDynamic()==false){
			m_permaSolidParticles.push_back(p);
		}else if(m_geoms[solidGeomID].m_geom->IsDynamic()==true){
			m_solidParticles.push_back(p);
		}else{
			delete p;
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

