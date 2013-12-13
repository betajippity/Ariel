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

}

scene::~scene(){

}

void scene::addSolidObject(objCore::objContainer* object){
	solidObjects.push_back(object);
	openvdb::FloatGrid::Ptr objectSDF;
	meshToSDF(objectSDF, object);
	if(solidObjects.size()==1){
		solidSDF = objectSDF->deepCopy();
	}else{
		openvdb::tools::csgUnion(*solidSDF, *objectSDF);
	}
	objectSDF->clear();
}

void scene::addLiquidObject(objCore::objContainer* object){
	liquidObjects.push_back(object);
	openvdb::FloatGrid::Ptr objectSDF;
	meshToSDF(objectSDF, object);
	if(liquidObjects.size()==1){
		liquidSDF = objectSDF->deepCopy();
	}else{
		openvdb::tools::csgUnion(*liquidSDF, *objectSDF);
	}
	objectSDF->clear();
}

vector<objCore::objContainer*>& scene::getSolidObjects(){
	return solidObjects;
}

vector<objCore::objContainer*>& scene::getLiquidObjects(){
	return liquidObjects;
}

void scene::meshToSDF(openvdb::FloatGrid::Ptr& grid, objCore::objContainer* mesh){
	openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(.25f);
	//copy vertices into vdb format
	vector<openvdb::Vec3s> vdbpoints;
	for(int i=0; i<mesh->getObj()->numberOfVertices; i++){
		vec3 vertex = mesh->getObj()->vertices[i];
		openvdb::Vec3s vdbvertex(vertex.x, vertex.y, vertex.z);
		vdbpoints.push_back(transform->worldToIndex(vdbvertex));
	}
	//copy faces into vdb format
	vector<openvdb::Vec4I> vdbpolys;
	for(int i=0; i<mesh->getObj()->numberOfPolys; i++){
		vec4 poly = mesh->getObj()->polyVertexIndices[i];
		openvdb::Vec4I vdbpoly((int)poly[0]-1, (int)poly[1]-1, (int)poly[2]-1, (int)poly[3]-1);
		if((int)poly[0]==(int)poly[3] || (int)poly[3]<0){
			vdbpoly[3] = openvdb::util::INVALID_IDX;
		}
		vdbpolys.push_back(vdbpoly);
	}	
	 //call vdb tools for creating level set
	openvdb::tools::MeshToVolume<openvdb::FloatGrid> sdfmaker(transform);
	sdfmaker.convertToLevelSet(vdbpoints, vdbpolys);
	grid = sdfmaker.distGridPtr();
	//cleanup
	vdbpoints.clear();
	vdbpolys.clear();
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

	fluidCore::floatgrid* solidGrid = new fluidCore::floatgrid(solidSDF);
	fluidCore::floatgrid* liquidGrid = new fluidCore::floatgrid(liquidSDF);

	if(type==FLUID){
		vec3 worldpos = pos*scale;
		if(liquidGrid->getInterpolatedCell(worldpos)<thickness){
			inside = true;
		}
		//if particles are in a wall, don't generate them
		if(solidGrid->getInterpolatedCell(worldpos)<thickness){ 
			inside = false; 
		}	

	}else if(type==SOLID){
		vec3 worldpos = pos*scale;
		if(solidGrid->getInterpolatedCell(worldpos)<thickness){
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
