// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: levelset.cpp
// Implements levelset.hpp

#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/tools/VolumeToSpheres.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/util/NullInterrupter.h>
#include "levelset.hpp"

namespace fluidCore{

levelset::levelset(){
	openvdb::initialize();
	vdbgrid = openvdb::FloatGrid::create(0.0f);
}

levelset::~levelset(){
	vdbgrid->clear();
	vdbgrid.reset();
}

levelset::levelset(objCore::objContainer* mesh){
	openvdb::math::Transform::Ptr transform=openvdb::math::Transform::createLinearTransform(.25f);
	//copy vertices into vdb format
	std::vector<openvdb::Vec3s> vdbpoints;
	for(int i=0; i<mesh->getObj()->numberOfVertices; i++){
		glm::vec3 vertex = mesh->getObj()->vertices[i];
		openvdb::Vec3s vdbvertex(vertex.x, vertex.y, vertex.z);
		vdbpoints.push_back(transform->worldToIndex(vdbvertex));
	}
	//copy faces into vdb format
	std::vector<openvdb::Vec4I> vdbpolys;
	for(int i=0; i<mesh->getObj()->numberOfPolys; i++){
		glm::vec4 poly = mesh->getObj()->polyVertexIndices[i];
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

levelset::levelset(std::vector<particle*>& particles, float maxdimension){
	vdbgrid = openvdb::createLevelSet<openvdb::FloatGrid>();
	openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*vdbgrid);
	raster.setGrainSize(1);
	raster.setRmin(.01f);

	particleList plist(particles, maxdimension);
	// raster.rasterizeSpheres(plist);
	raster.rasterizeTrails(plist);
	raster.finalize();
}

void levelset::writeObjToFile(std::string filename){
	openvdb::tools::VolumeToMesh vdbmesher(0,.05f);
	vdbmesher(*getVDBGrid());

	//get all mesh points and dump to glm::vec3 std::vector
	int vdbPointsCount = vdbmesher.pointListSize();
	openvdb::tools::PointList& vdbPoints = vdbmesher.pointList();
	
	std::vector<glm::vec3> points;
	points.reserve(vdbPointsCount);
	for(int i=0; i<vdbPointsCount; i++){
		glm::vec3 v(vdbPoints[i][0], vdbPoints[i][1], vdbPoints[i][2]);
		points.push_back(v);
	}

	//get all mesh faces and dump to glm::vec4 std::vector
	int vdbFacesCount = vdbmesher.polygonPoolListSize();
	openvdb::tools::PolygonPoolList& vdbFaces = vdbmesher.polygonPoolList();

	int facesCount = 0;
	for(int i=0; i<vdbFacesCount; i++){
		facesCount = facesCount + vdbFaces[i].numQuads() + vdbFaces[i].numTriangles();
	}

	std::vector<glm::vec4> faces;
	faces.reserve(facesCount);

	int test = 5;

	for(int i=0; i<vdbFacesCount; i++){
		int count = vdbFaces[i].numQuads();
		for(int j=0; j<count; j++){
			openvdb::Vec4I vdbface = vdbFaces[i].quad(j);
			glm::vec4 f(vdbface[0]+1, vdbface[1]+1, vdbface[2]+1, vdbface[3]+1);
			faces.push_back(f);
		}
		count = vdbFaces[i].numTriangles();
		for(int j=0; j<count; j++){
			openvdb::Vec3I vdbface = vdbFaces[i].triangle(j);
			glm::vec4 f(vdbface[0]+1, vdbface[1]+1, vdbface[2]+1, -1);
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

void levelset::projectPointsToSurface(std::vector<glm::vec3>& points){
	std::vector<openvdb::Vec3R> vdbpoints;
	std::vector<float> distances;
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
		points[i] = glm::vec3(vdbpoints[i][0], vdbpoints[i][1], vdbpoints[i][2]);
	}
}

void levelset::writeVDBGridToFile(std::string filename){
	openvdb::io::File file(filename);
	openvdb::GridPtrVec grids;
    grids.push_back(vdbgrid);
    file.write(grids);
    file.close();
}

float levelset::getInterpolatedCell(const glm::vec3& index){
	return getInterpolatedCell(index.x, index.y, index.z);
}

float levelset::getInterpolatedCell(const float& x, const float& y, const float& z){
	float value;
	GetInterpolatedCellLock.lock();
	{
		openvdb::Vec3f p(x,y,z);
		openvdb::tools::GridSampler<openvdb::FloatTree, openvdb::tools::BoxSampler> interpolator(
													  vdbgrid->constTree(), vdbgrid->transform());
		value = interpolator.wsSample(p);
	}
	GetInterpolatedCellLock.unlock();
	return value;
}

float levelset::getCell(const glm::vec3& index){
	return getCell((int)index.x, (int)index.y, (int)index.z);
}

float levelset::getCell(const int& x, const int& y, const int& z){
	openvdb::Coord coord = openvdb::Coord(x,y,z);
	openvdb::FloatGrid::Accessor accessor = vdbgrid->getAccessor();
	float cell = accessor.getValue(coord);
	return cell;
}

void levelset::setCell(const glm::vec3& index, const float& value){
	setCell((int)index.x, (int)index.y, (int)index.z, value);
}

void levelset::setCell(const int& x, const int& y, const int& z, const float& value){
	SetCellLock.lock();
	{
		openvdb::Coord coord = openvdb::Coord(x,y,z);
		openvdb::FloatGrid::Accessor accessor = vdbgrid->getAccessor();
		accessor.setValue(coord, value);
	}
	SetCellLock.unlock();
}

openvdb::FloatGrid::Ptr& levelset::getVDBGrid(){
	return vdbgrid;
}

void levelset::merge(levelset& ls){
	openvdb::FloatGrid::Ptr objectSDF = ls.getVDBGrid()->deepCopy();
	openvdb::tools::csgUnion(*vdbgrid, *objectSDF);
	objectSDF->clear();
}
	
void levelset::copy(levelset& ls){
	vdbgrid = ls.getVDBGrid()->deepCopy();
}
}
