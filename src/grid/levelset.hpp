// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: levelset.hpp
// Class that extends floatgrid for representing level sets

#ifndef LEVELSET_HPP
#define LEVELSET_HPP

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Composite.h>
#include "macgrid.inl"
#include "../geom/geomlist.hpp"

namespace fluidCore {
//====================================
// Class Declarations
//====================================

class ParticleList{ //used for VDB particle to level set construction
	public:
		ParticleList(){ }

		ParticleList(std::vector<Particle*> plist, float maxdimension){
			m_particles = plist;
			m_maxdimension = maxdimension;
		}

		~ParticleList(){ }

		int size() const { 
			return m_particles.size(); 
		}

		void getPos(size_t n, openvdb::Vec3R& pos) const {
			pos = openvdb::Vec3f(m_particles[n]->m_p.x*m_maxdimension, 
								 m_particles[n]->m_p.y*m_maxdimension, 
								 m_particles[n]->m_p.z*m_maxdimension);
		}

		void getPosRad(size_t n, openvdb::Vec3R& pos, openvdb::Real& rad) const {
		    pos = openvdb::Vec3f(m_particles[n]->m_p.x*m_maxdimension, 
		    					 m_particles[n]->m_p.y*m_maxdimension, 
		    					 m_particles[n]->m_p.z*m_maxdimension);
		    rad = m_particles[n]->m_density;
		    rad = .5f;
		    if(m_particles[n]->m_invalid){
		    	rad = 0.0f;
		    }
		}

		void getPosRadVel(size_t n, openvdb::Vec3R& pos, openvdb::Real& rad, 
						  openvdb::Vec3R& vel) const {
			pos = openvdb::Vec3f(m_particles[n]->m_p.x*m_maxdimension, 
								 m_particles[n]->m_p.y*m_maxdimension, 
								 m_particles[n]->m_p.z*m_maxdimension);
		    rad = m_particles[n]->m_density;
		    rad = .5f;
		    vel = openvdb::Vec3f(m_particles[n]->m_u.x, m_particles[n]->m_u.y, 
		    					 m_particles[n]->m_u.z);
		    if(m_particles[n]->m_invalid){
		    	rad = 0.0f;
		    }
	    }

		void getAtt(size_t n, openvdb::Index32& att) const { att = n; }
	private:
		std::vector<Particle*>		m_particles;
		float						m_maxdimension;
};

class LevelSet{
	public:
		//Initializers
		LevelSet();
		// levelset(openvdb::FloatGrid::Ptr grid);
		LevelSet(objCore::Obj* mesh);
		LevelSet(objCore::Obj* mesh, const glm::mat4& m);
		LevelSet(objCore::InterpolatedObj* animmesh, const float& interpolation, 
				 const glm::mat4& m);
		LevelSet(std::vector<Particle*>& particles, float maxdimension);
		~LevelSet();

		//Cell accessors and setters and whatever
		float GetCell(const glm::vec3& index);
		float GetCell(const int& x, const int& y, const int& z);

		void SetCell(const glm::vec3& index, const float& value);
		void SetCell(const int& x, const int& y, const int& z, const float& value);

		float GetInterpolatedCell(const glm::vec3& index);
		float GetInterpolatedCell(const float& x, const float& y, const float& z);

		openvdb::FloatGrid::Ptr& GetVDBGrid();

		void Merge(LevelSet& ls);
		void Copy(LevelSet& ls);

		void ProjectPointsToSurface(std::vector<Particle*>& particles, const float& pscale);

		void WriteObjToFile(std::string filename);
		void WriteVDBGridToFile(std::string filename);

	protected:
		void LevelSetFromAnimMesh(objCore::InterpolatedObj* animmesh, const float& interpolation, 
								  const glm::mat4& m);
		void LevelSetFromMesh(objCore::Obj* mesh, const glm::mat4& m);

		openvdb::FloatGrid::Ptr		m_vdbgrid;

		tbb::mutex					m_getInterpolatedCellLock;
		tbb::mutex					m_setCellLock;
};
}

#endif
