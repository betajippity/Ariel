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

LevelSet::LevelSet(){
	openvdb::initialize();
	m_vdbgrid = openvdb::FloatGrid::create(0.0f);
}

LevelSet::~LevelSet(){
	m_vdbgrid->clear();
	m_vdbgrid.reset();
}

LevelSet::LevelSet(objCore::Obj* mesh){
	openvdb::math::Transform::Ptr transform=openvdb::math::Transform::createLinearTransform(.25f);
	//copy vertices into vdb format
	std::vector<openvdb::Vec3s> vdbpoints;
	for(unsigned int i=0; i<mesh->m_numberOfVertices; i++){
		glm::vec3 vertex = mesh->m_vertices[i];
		openvdb::Vec3s vdbvertex(vertex.x, vertex.y, vertex.z);
		vdbpoints.push_back(transform->worldToIndex(vdbvertex));
	}
	//copy faces into vdb format
	std::vector<openvdb::Vec4I> vdbpolys;
	for(unsigned int i=0; i<mesh->m_numberOfPolys; i++){
		glm::uvec4 poly = mesh->m_polyVertexIndices[i];
		openvdb::Vec4I vdbpoly(poly[0]-1, poly[1]-1, poly[2]-1, poly[3]-1);
		if(poly[0]==poly[3] || poly[3]<0){
			vdbpoly[3] = openvdb::util::INVALID_IDX;
		}
		vdbpolys.push_back(vdbpoly);
	}	
	 //call vdb tools for creating level set
	openvdb::tools::MeshToVolume<openvdb::FloatGrid> sdfmaker(transform);
	sdfmaker.convertToLevelSet(vdbpoints, vdbpolys);
	m_vdbgrid = sdfmaker.distGridPtr();
	//cleanup
	vdbpoints.clear();
	vdbpolys.clear();
}

LevelSet::LevelSet(std::vector<Particle*>& particles, float maxdimension){
	m_vdbgrid = openvdb::createLevelSet<openvdb::FloatGrid>();
	openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*m_vdbgrid);
	raster.setGrainSize(1);
	raster.setRmin(.01f);

	ParticleList plist(particles, maxdimension);
	// raster.rasterizeSpheres(plist);
	raster.rasterizeTrails(plist);
	raster.finalize();
}

void LevelSet::WriteObjToFile(std::string filename){
	openvdb::tools::VolumeToMesh vdbmesher(0,.05f);
	vdbmesher(*GetVDBGrid());

	//get all mesh points and dump to glm::vec3 std::vector
	unsigned int vdbPointsCount = vdbmesher.pointListSize();
	openvdb::tools::PointList& vdbPoints = vdbmesher.pointList();
	
	std::vector<glm::vec3> points;
	points.reserve(vdbPointsCount);
	for(unsigned int i=0; i<vdbPointsCount; i++){
		glm::vec3 v(vdbPoints[i][0], vdbPoints[i][1], vdbPoints[i][2]);
		points.push_back(v);
	}

	//get all mesh faces and dump to glm::vec4 std::vector
	unsigned int vdbFacesCount = vdbmesher.polygonPoolListSize();
	openvdb::tools::PolygonPoolList& vdbFaces = vdbmesher.polygonPoolList();

	unsigned int facesCount = 0;
	for(unsigned int i=0; i<vdbFacesCount; i++){
		facesCount = facesCount + vdbFaces[i].numQuads() + vdbFaces[i].numTriangles();
	}

	std::vector<glm::uvec4> faces;
	faces.reserve(facesCount);

	unsigned int test = 5;

	for(unsigned int i=0; i<vdbFacesCount; i++){
		unsigned int count = vdbFaces[i].numQuads();
		for(unsigned int j=0; j<count; j++){
			openvdb::Vec4I vdbface = vdbFaces[i].quad(j);
			glm::uvec4 f(vdbface[0]+1, vdbface[1]+1, vdbface[2]+1, vdbface[3]+1);
			faces.push_back(f);
		}
		count = vdbFaces[i].numTriangles();
		for(unsigned int j=0; j<count; j++){
			openvdb::Vec3I vdbface = vdbFaces[i].triangle(j);
			glm::uvec4 f(vdbface[0]+1, vdbface[1]+1, vdbface[2]+1, -1);
			faces.push_back(f);
		}
	}

	//pack points and faces into objcontainer and write
	objCore::Obj* mesh = new objCore::Obj();
    mesh->m_numberOfVertices = points.size();
    mesh->m_vertices = &points[0];
    mesh->m_numberOfNormals = 0;
    mesh->m_normals = NULL;
    mesh->m_numberOfUVs = 0;
    mesh->m_uvs = NULL;
    mesh->m_numberOfPolys = faces.size();
    mesh->m_polyVertexIndices = &faces[0];
    mesh->m_polyNormalIndices = NULL;
    mesh->m_polyUVIndices = NULL;
    
	mesh->m_keep = true;

	mesh->WriteObj(filename);

	delete mesh;
}

void LevelSet::ProjectPointsToSurface(std::vector<glm::vec3>& points){
	std::vector<openvdb::Vec3R> vdbpoints;
	std::vector<float> distances;
	unsigned int pointsCount = points.size();
	vdbpoints.reserve(pointsCount);
	distances.reserve(pointsCount);
	for(unsigned int i=0; i<pointsCount; i++){
		openvdb::Vec3s vdbvertex(points[i].x, points[i].y, points[i].z);
		vdbpoints.push_back(vdbvertex);
	}
	openvdb::tools::ClosestSurfacePoint<openvdb::FloatGrid> csp;
	openvdb::util::NullInterrupter n;
	csp.initialize(*m_vdbgrid, 0.0f, &n);
	csp.searchAndReplace(vdbpoints, distances);
	for(unsigned int i=0; i<pointsCount; i++){
		points[i] = glm::vec3(vdbpoints[i][0], vdbpoints[i][1], vdbpoints[i][2]);
	}
}

void LevelSet::WriteVDBGridToFile(std::string filename){
	openvdb::io::File file(filename);
	openvdb::GridPtrVec grids;
    grids.push_back(m_vdbgrid);
    file.write(grids);
    file.close();
}

float LevelSet::GetInterpolatedCell(const glm::vec3& index){
	return GetInterpolatedCell(index.x, index.y, index.z);
}

float LevelSet::GetInterpolatedCell(const float& x, const float& y, const float& z){
	float value;
	m_getInterpolatedCellLock.lock();
	{
		openvdb::Vec3f p(x,y,z);
		openvdb::tools::GridSampler<openvdb::FloatTree, openvdb::tools::BoxSampler> interpolator(
													  					m_vdbgrid->constTree(), 
													  					m_vdbgrid->transform());
		value = interpolator.wsSample(p);
	}
	m_getInterpolatedCellLock.unlock();
	return value;
}

float LevelSet::GetCell(const glm::vec3& index){
	return GetCell((int)index.x, (int)index.y, (int)index.z);
}

float LevelSet::GetCell(const int& x, const int& y, const int& z){
	openvdb::Coord coord = openvdb::Coord(x,y,z);
	openvdb::FloatGrid::Accessor accessor = m_vdbgrid->getAccessor();
	float cell = accessor.getValue(coord);
	return cell;
}

void LevelSet::SetCell(const glm::vec3& index, const float& value){
	SetCell((int)index.x, (int)index.y, (int)index.z, value);
}

void LevelSet::SetCell(const int& x, const int& y, const int& z, const float& value){
	m_setCellLock.lock();
	{
		openvdb::Coord coord = openvdb::Coord(x,y,z);
		openvdb::FloatGrid::Accessor accessor = m_vdbgrid->getAccessor();
		accessor.setValue(coord, value);
	}
	m_setCellLock.unlock();
}

openvdb::FloatGrid::Ptr& LevelSet::GetVDBGrid(){
	return m_vdbgrid;
}

void LevelSet::Merge(LevelSet& ls){
	openvdb::FloatGrid::Ptr objectSDF = ls.GetVDBGrid()->deepCopy();
	openvdb::tools::csgUnion(*m_vdbgrid, *objectSDF);
	objectSDF->clear();
}
	
void LevelSet::Copy(LevelSet& ls){
	m_vdbgrid = ls.GetVDBGrid()->deepCopy();
}
}
