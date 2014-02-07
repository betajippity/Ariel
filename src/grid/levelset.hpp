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

		particleList(vector<particle*> plist, float maxdimension){
			particles = plist;
			this->maxdimension = maxdimension;
		}

		~particleList(){ }

		int size() const { 
			return particles.size(); 
		}

		void getPos(size_t n, openvdb::Vec3R& pos) const {
			pos = openvdb::Vec3f(particles[n]->p.x*maxdimension, particles[n]->p.y*maxdimension, 
								 particles[n]->p.z*maxdimension);
		}

		void getPosRad(size_t n, openvdb::Vec3R& pos, openvdb::Real& rad) const {
		    pos = openvdb::Vec3f(particles[n]->p.x*maxdimension, particles[n]->p.y*maxdimension, 
		    					 particles[n]->p.z*maxdimension);
		    rad = particles[n]->density;
		    rad = .5f;
		    if(particles[n]->invalid){
		    	rad = 0.0f;
		    }
		}

		void getPosRadVel(size_t n, openvdb::Vec3R& pos, openvdb::Real& rad, openvdb::Vec3R& vel) const {
			pos = openvdb::Vec3f(particles[n]->p.x*maxdimension, particles[n]->p.y*maxdimension, 
								 particles[n]->p.z*maxdimension);
		    rad = particles[n]->density;
		    rad = .5f;
		    vel = openvdb::Vec3f(particles[n]->u.x, particles[n]->u.y, particles[n]->u.z);
		    if(particles[n]->invalid){
		    	rad = 0.0f;
		    }
	    }

		void getAtt(size_t n, openvdb::Index32& att) const { att = n; }
	private:
		vector<particle*> particles;
		float maxdimension;
};

class levelset: public floatgrid{
	public:
		//Initializers
		levelset();
		// levelset(openvdb::FloatGrid::Ptr grid);
		levelset(objCore::objContainer* mesh);
		levelset(vector<particle*>& particles, float maxdimension);
		~levelset();

		void merge(levelset& ls);
		void copy(levelset& ls);

		void projectPointsToSurface(vector<vec3>& points);

		void writeObjToFile(string filename);
};
}

#endif
