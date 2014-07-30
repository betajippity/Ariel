// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: aabb.hpp. Adapted from Takua Render.
// Axis aligned bounding box stuff

#ifndef AABB_HPP
#define AABB_HPP

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#define SHARED __shared__
#else
#define HOST
#define DEVICE
#define SHARED
#endif

#include "../ray/ray.hpp"
#include "../utilities/utilities.h"

namespace spaceCore {

//====================================
// Class Declarations
//====================================

class Aabb {
    public:
        Aabb();
        Aabb(const glm::vec3& min, const glm::vec3& max, const int& id);
        Aabb(const glm::vec3& min, const glm::vec3& max, const glm::vec3& centroid, const int& id);
        ~Aabb();

        void SetContents(const glm::vec3& min, const glm::vec3& max, const glm::vec3& centroid,
                         const int& id);
        void ExpandAabb(const glm::vec3& exMin, const glm::vec3& exMax);
        double CalculateSurfaceArea();
        //"fast" as in only returns distance and not a full intersection           
        HOST DEVICE float FastIntersectionTest(const rayCore::Ray& r); 

        glm::vec3   m_min;
        glm::vec3   m_max;
        glm::vec3   m_centroid;
        int         m_id; //what this id actually means is entirely dependent on the context
};
}

#endif
