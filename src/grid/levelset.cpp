// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: levelset.cpp
// Implements levelset.hpp

#include "levelset.hpp"
#include <openvdb/tools/ParticlesToLevelSet.h>

using namespace fluidCore;
using namespace utilityCore;

levelset::levelset(){
	type = VDB;
}

levelset::~levelset(){

}

levelset::levelset(objCore::objContainer* mesh){
	type = VDB;
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
	vdbgrid = sdfmaker.distGridPtr();
	//cleanup
	vdbpoints.clear();
	vdbpolys.clear();
}

levelset::levelset(vector<particle*>& particles){
	type = VDB;
	vdbgrid = openvdb::createLevelSet<openvdb::FloatGrid>();
	openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*vdbgrid);
	raster.setGrainSize(1);
	raster.setRmin(.01f);

	particleList plist(particles);
	raster.rasterizeSpheres(plist);
	raster.finalize();
}

void levelset::merge(levelset& ls){
	openvdb::FloatGrid::Ptr objectSDF = ls.getVDBGrid()->deepCopy();
	openvdb::tools::csgUnion(*vdbgrid, *objectSDF);
	objectSDF->clear();
}
	
void levelset::copy(levelset& ls){
	vdbgrid = ls.getVDBGrid()->deepCopy();
}
