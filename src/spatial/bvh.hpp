// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: bvh.hpp. Adapted from Takua Render.
// Spatial BVH acceleration structure

#ifndef BVH_HPP
#define BVH_HPP

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#define SHARED __shared__
#else
#define HOST
#define DEVICE
#define SHARED
#endif

#include "aabb.hpp"
#include "spatial.hpp"
#include "../ray/ray.hpp"
#include "../utilities/utilities.h"

enum Axis{axis_x, axis_y, axis_z};

namespace spaceCore {

//====================================
// Struct Declarations
//====================================

struct BvhNode {
    Aabb m_bounds;
    unsigned int m_referenceOffset; //offset into bvh-wide reference index 
    unsigned int m_numberOfReferences;
    unsigned int m_left;
    unsigned int m_right;
    unsigned int m_nodeid;

    BvhNode(){
        m_bounds = Aabb();
        m_referenceOffset = 0;
        m_numberOfReferences = 0;
        m_left = 0;
        m_right = 0;
    }

    HOST DEVICE bool IsLeaf(){
        if(m_left==0 && m_right==0){
            return true;
        }else{
            return false;
        }
    }

    HOST DEVICE float FastIntersectionTest(const rayCore::Ray& r){
        return m_bounds.FastIntersectionTest(r);
    }
};

//====================================
// Class Declarations
//===================================

template <typename T> class Bvh {
    public:
        Bvh(T basegeom);
        Bvh();
        ~Bvh();

        void BuildBvh(const unsigned int& maxDepth);
        HOST DEVICE void Traverse(const rayCore::Ray& r, TraverseAccumulator& result);

        BvhNode*                    m_nodes;
        unsigned int                m_numberOfNodes;
        unsigned int*               m_referenceIndices;
        unsigned int                m_numberOfReferenceIndices;
        unsigned int                m_id;
        unsigned int                m_depth;
        T                           m_basegeom;

    private:
        float FindBestSplit(BvhNode& Node, std::vector<unsigned int>& references,
                            const Axis& direction, const unsigned int& quality, Aabb* aabbs,
                            unsigned int& leftCount, unsigned int& rightCount);
        Axis FindLongestAxis(const Aabb& box);
        double CalculateSplitCost(const float& split, const BvhNode& node, 
                                  std::vector<unsigned int>& references, const Axis& direction,
                                  Aabb* aabbs, unsigned int& leftCount, unsigned int& rightCount);
};
}

#include "bvh.inl"

#endif
