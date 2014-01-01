// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: levelset.hpp
// Class that extends floatgrid for representing level sets

#ifndef LEVELSET_HPP
#define LEVELSET_HPP

#include "floatgrid.hpp"
#include "macgrid.inl"
#include "../geom/geom.inl"
#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Composite.h>

using namespace std;
using namespace glm;

namespace fluidCore {
//====================================
// Class Declarations
//====================================

class particleList{ //used for VDB particle to level set construction
	public:
		particleList(){ }

		particleList(vector<particle*> plist){
			particles = plist;
		}

		~particleList(){ }

		int size() const { 
			return particles.size(); 
		}

		void getPos(size_t n, openvdb::Vec3R& pos) const {
			pos = openvdb::Vec3f(particles[n]->p.x*32.0f, particles[n]->p.y*32.0f, particles[n]->p.z*32.0f);
		}

		void getPosRad(size_t n, openvdb::Vec3R& pos, openvdb::Real& rad) const {
		    pos = openvdb::Vec3f(particles[n]->p.x*32.0f, particles[n]->p.y*32.0f, particles[n]->p.z*32.0f);
		    rad = particles[n]->density;
		}

		void getPosRadVel(size_t n, openvdb::Vec3R& pos, openvdb::Real& rad, openvdb::Vec3R& vel) const {
			pos = openvdb::Vec3f(particles[n]->p.x*32.0f, particles[n]->p.y*32.0f, particles[n]->p.z*32.0f);
		    rad = particles[n]->density;
		    vel = openvdb::Vec3f(particles[n]->u.x, particles[n]->u.y, particles[n]->u.z);
	    }

		void getAtt(size_t n, openvdb::Index32& att) const { att = n; }
	private:
		vector<particle*> particles;
};

class levelset: public floatgrid{
	public:
		//Initializers
		levelset();
		// levelset(openvdb::FloatGrid::Ptr grid);
		levelset(objCore::objContainer* mesh);
		levelset(vector<particle*>& particles);
		~levelset();

		void merge(levelset& ls);
		void copy(levelset& ls);
};
}

#endif
