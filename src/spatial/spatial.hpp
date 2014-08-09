// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: spatial.hpp. Adapted from Takua Render.
// Spatial acceleration stuff

#ifndef SPATIAL_HPP
#define SPATIAL_HPP

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
//===================================

class TraverseAccumulator {
    public:
        HOST DEVICE TraverseAccumulator();
        HOST DEVICE TraverseAccumulator(const glm::vec3& origin);
        HOST DEVICE ~TraverseAccumulator();
        
        HOST DEVICE virtual void RecordIntersection(const rayCore::Intersection& intersect, 
                                                    const unsigned int& nodeid);
        HOST DEVICE virtual void Transform(const glm::mat4& m);

        rayCore::Intersection   m_intersection;
        unsigned int            m_nodeid;
        glm::vec3               m_origin;
};

class HitCountTraverseAccumulator: public TraverseAccumulator {
    public:
        HOST DEVICE HitCountTraverseAccumulator();
        HOST DEVICE HitCountTraverseAccumulator(const glm::vec3& origin);
        HOST DEVICE ~HitCountTraverseAccumulator();

        HOST DEVICE void RecordIntersection(const rayCore::Intersection& intersect, 
                                            const unsigned int& nodeid);
        HOST DEVICE void Transform(const glm::mat4& m);

        rayCore::Intersection   m_intersection;
        unsigned int            m_nodeid;  
        unsigned int            m_numberOfHits;
        glm::vec3               m_origin;
};

class DebugTraverseAccumulator: public TraverseAccumulator {
    public:
        DebugTraverseAccumulator();
        ~DebugTraverseAccumulator();

        void RecordIntersection(const rayCore::Intersection& intersect,
                                const unsigned int& nodeid);
        void Transform(const glm::mat4& m);

        std::vector<rayCore::Intersection> m_intersections;
        std::vector<unsigned int> m_nodeids;
};
}

#endif
