// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: mesh.hpp
// Mesh class inherits geomContainer for obj meshes

#ifndef MESH_HPP
#define MESH_HPP

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#define SHARED __shared__
#else
#define HOST
#define DEVICE
#define SHARED
#endif

#include "geom.hpp"
#include "obj/obj.hpp"
#include "../spatial/bvh.hpp"
#include "../utilities/utilities.h"

namespace geomCore {

//====================================
// Struct Declarations
//====================================

//Only used as a temp data transfer container for CUDA memory setup
struct AnimatedMeshContainerData {
    spaceCore::Bvh<objCore::InterpolatedObj>**      m_meshFrames;
    GeomTransform**                                 m_geomTransforms;
    unsigned int                                    m_numberOfFrames;
    unsigned int                                    m_id; 
    unsigned int                                    m_frameOffset;
    unsigned int                                    m_frameInterval;
};

struct MeshContainerData {
    spaceCore::Bvh<objCore::Obj>**                  m_meshFrames;
    GeomTransform**                                 m_geomTransforms;
    unsigned int                                    m_numberOfFrames;
    unsigned int                                    m_id; 
    unsigned int                                    m_frameOffset;
    unsigned int                                    m_frameInterval;
};

//====================================
// Class Declarations
//====================================

class MeshContainer: public virtual GeomInterface {
    public:
        HOST DEVICE MeshContainer();
        HOST DEVICE MeshContainer(const unsigned int& numberOfFrames,
                                  const unsigned int& frameOffset,
                                  const unsigned int& frameInterval,
                                  const bool& prePersist,
                                  const bool& postPersist,
                                  GeomTransform** geomTransforms,
                                  spaceCore::Bvh<objCore::Obj>** meshFrames);
        HOST DEVICE MeshContainer(MeshContainerData data); 
        HOST DEVICE ~MeshContainer();

        HOST DEVICE GeomType GetType();
        HOST DEVICE unsigned int GetID();
        HOST DEVICE void Intersect(const rayCore::Ray& r, spaceCore::TraverseAccumulator& result);

        spaceCore::Bvh<objCore::Obj>**                  m_meshFrames;
        GeomTransform**                                 m_geomTransforms;
        unsigned int                                    m_numberOfFrames;
        unsigned int                                    m_id;
        unsigned int                                    m_frameOffset;
        unsigned int                                    m_frameInterval;
        bool                                            m_prePersist;
        bool                                            m_postPersist;
};

class AnimatedMeshContainer: public virtual GeomInterface {
    public:
        HOST DEVICE AnimatedMeshContainer();
        HOST DEVICE AnimatedMeshContainer(const unsigned int& numberOfFrames,
                                          const unsigned int& frameOffset,
                                          const unsigned int& frameInterval,
                                          const bool& prePersist,
                                          const bool& postPersist,
                                          GeomTransform** geomTransforms,
                                          spaceCore::Bvh<objCore::InterpolatedObj>** meshFrames);
        HOST DEVICE AnimatedMeshContainer(AnimatedMeshContainerData data); 
        HOST DEVICE ~AnimatedMeshContainer();

        HOST DEVICE GeomType GetType();
        HOST DEVICE unsigned int GetID();
        HOST DEVICE void Intersect(const rayCore::Ray& r, spaceCore::TraverseAccumulator& result);

        spaceCore::Bvh<objCore::InterpolatedObj>**      m_meshFrames;
        GeomTransform**                                 m_geomTransforms;
        unsigned int                                    m_numberOfFrames;
        unsigned int                                    m_id;
        unsigned int                                    m_frameOffset;
        unsigned int                                    m_frameInterval;
        bool                                            m_prePersist;
        bool                                            m_postPersist;
};
}

#endif
