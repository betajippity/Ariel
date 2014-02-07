// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: levelset.cpp
// Implements levelset.hpp

#include "levelset.hpp"
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/tools/VolumeToSpheres.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/util/NullInterrupter.h>

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

levelset::levelset(vector<particle*>& particles, float maxdimension){
	type = VDB;
	vdbgrid = openvdb::createLevelSet<openvdb::FloatGrid>();
	openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*vdbgrid);
	raster.setGrainSize(1);
	raster.setRmin(.01f);

	particleList plist(particles, maxdimension);
	// raster.rasterizeSpheres(plist);
	raster.rasterizeTrails(plist);
	raster.finalize();
}

void levelset::writeObjToFile(string filename){
	openvdb::tools::VolumeToMesh vdbmesher(0,.05f);
	vdbmesher(*getVDBGrid());

	//get all mesh points and dump to vec3 vector
	int vdbPointsCount = vdbmesher.pointListSize();
	openvdb::tools::PointList& vdbPoints = vdbmesher.pointList();
	
	vector<vec3> points;
	points.reserve(vdbPointsCount);
	for(int i=0; i<vdbPointsCount; i++){
		vec3 v(vdbPoints[i][0], vdbPoints[i][1], vdbPoints[i][2]);
		points.push_back(v);
	}

	//get all mesh faces and dump to vec4 vector
	int vdbFacesCount = vdbmesher.polygonPoolListSize();
	openvdb::tools::PolygonPoolList& vdbFaces = vdbmesher.polygonPoolList();

	int facesCount = 0;
	for(int i=0; i<vdbFacesCount; i++){
		facesCount = facesCount + vdbFaces[i].numQuads() + vdbFaces[i].numTriangles();
	}

	vector<vec4> faces;
	faces.reserve(facesCount);

	int test = 5;

	for(int i=0; i<vdbFacesCount; i++){
		int count = vdbFaces[i].numQuads();
		for(int j=0; j<count; j++){
			openvdb::Vec4I vdbface = vdbFaces[i].quad(j);
			vec4 f(vdbface[0]+1, vdbface[1]+1, vdbface[2]+1, vdbface[3]+1);
			faces.push_back(f);
		}
		count = vdbFaces[i].numTriangles();
		for(int j=0; j<count; j++){
			openvdb::Vec3I vdbface = vdbFaces[i].triangle(j);
			vec4 f(vdbface[0]+1, vdbface[1]+1, vdbface[2]+1, -1);
			faces.push_back(f);
		}
	}

	//pack points and faces into objcontainer and write
	objCore::obj* mesh = objCore::createObj(points.size(), &points[0], 0, NULL, 0, NULL, 
											faces.size(), &faces[0], NULL, NULL);
	objCore::objContainer* meshContainer = new objCore::objContainer(mesh);
	meshContainer->keepObj(true);

	meshContainer->writeObj(filename);

	delete mesh;
}

void levelset::projectPointsToSurface(vector<vec3>& points){
	vector<openvdb::Vec3R> vdbpoints;
	vector<float> distances;
	int pointsCount = points.size();
	vdbpoints.reserve(pointsCount);
	distances.reserve(pointsCount);
	for(int i=0; i<pointsCount; i++){
		openvdb::Vec3s vdbvertex(points[i].x, points[i].y, points[i].z);
		vdbpoints.push_back(vdbvertex);
	}
	openvdb::tools::ClosestSurfacePoint<openvdb::FloatGrid> csp;
	openvdb::util::NullInterrupter n;
	csp.initialize(*vdbgrid, 0.0f, &n);
	csp.searchAndReplace(vdbpoints, distances);
	for(int i=0; i<pointsCount; i++){
		points[i] = vec3(vdbpoints[i][0], vdbpoints[i][1], vdbpoints[i][2]);
	}
}

void levelset::merge(levelset& ls){
	openvdb::FloatGrid::Ptr objectSDF = ls.getVDBGrid()->deepCopy();
	openvdb::tools::csgUnion(*vdbgrid, *objectSDF);
	objectSDF->clear();
}
	
void levelset::copy(levelset& ls){
	vdbgrid = ls.getVDBGrid()->deepCopy();
}
