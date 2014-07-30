// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: mesh.cpp
// Implements mesh.hpp

#include "mesh.hpp"

namespace geomCore {

//====================================
// MeshContainer Class
//====================================

HOST DEVICE MeshContainer::MeshContainer(){
    m_meshFrames = NULL;
    m_geomFrames = NULL;
    m_prePersist = false;
    m_postPersist = false;
    m_frameInterval = 1;
    m_frameOffset = 0;
}

HOST DEVICE MeshContainer::MeshContainer(const unsigned int& numberOfFrames,
                                         const unsigned int& frameOffset,
                                         const unsigned int& frameInterval,
                                         const bool& prePersist,
                                         const bool& postPersist,
                                         GeomFrame** geomFrames,
                                         spaceCore::Bvh<objCore::Obj>** meshFrames){
    m_numberOfFrames = numberOfFrames;
    m_geomFrames = geomFrames;
    m_meshFrames = meshFrames;
    m_frameOffset = frameOffset;
    m_frameInterval = frameInterval;
    m_prePersist = prePersist;
    m_postPersist = postPersist;
}

HOST DEVICE MeshContainer::MeshContainer(MeshContainerData data){
    m_numberOfFrames = data.m_numberOfFrames;
    m_geomFrames = data.m_geomFrames;
    m_meshFrames = data.m_meshFrames;
    m_id = data.m_id;
}

HOST DEVICE MeshContainer::~MeshContainer(){

}

HOST DEVICE GeomType MeshContainer::GetType(){
    return MESH;
}

HOST DEVICE unsigned int MeshContainer::GetID(){
    return m_id;
}

HOST DEVICE void MeshContainer::Intersect(const rayCore::Ray& r, 
                                          spaceCore::TraverseAccumulator& result){
    //translate frame into local frame space
    float clampedFrame = (r.m_frame - float(m_frameOffset)) / float(m_frameInterval);
    if((clampedFrame<0.0f && m_prePersist==false) || 
       (clampedFrame>=m_numberOfFrames && m_postPersist==false)){
        return;
    }
    clampedFrame = glm::clamp(clampedFrame, 0.0f, float(m_numberOfFrames-1));
    //ceil/floor to get frame indices
    unsigned int upperFrame = glm::ceil(clampedFrame);
    unsigned int lowerFrame = glm::floor(clampedFrame);
    float lerpWeight = clampedFrame - float(lowerFrame);
    //grab relevant transforms and LERP to get interpolated transform, apply inverse to ray
    glm::vec3 interpT = m_geomFrames[lowerFrame]->m_translation * (1.0f-lerpWeight) + 
                        m_geomFrames[upperFrame]->m_translation * lerpWeight;
    glm::vec3 interpR = m_geomFrames[lowerFrame]->m_rotation * (1.0f-lerpWeight) + 
                        m_geomFrames[upperFrame]->m_rotation * lerpWeight;
    glm::vec3 interpS = m_geomFrames[lowerFrame]->m_scale * (1.0f-lerpWeight) + 
                        m_geomFrames[upperFrame]->m_scale * lerpWeight;
    rayCore::Ray transformedR = r.Transform(utilityCore::buildInverseTransformationMatrix(interpT,
                                                                                          interpR,
                                                                                          interpS));
    //run intersection, transform result back into worldspace
    m_meshFrames[lowerFrame]->Traverse(transformedR, result);
    result.Transform(utilityCore::buildTransformationMatrix(interpT, interpR, interpS));
}

//====================================
// AnimatedMeshContainer Class
//====================================

HOST DEVICE AnimatedMeshContainer::AnimatedMeshContainer(){
    m_meshFrames = NULL;
    m_geomFrames = NULL;
    m_prePersist = false;
    m_postPersist = false;
    m_frameInterval = 1;
    m_frameOffset = 0;
}

HOST DEVICE AnimatedMeshContainer::AnimatedMeshContainer(const unsigned int& numberOfFrames,
                                                         const unsigned int& frameOffset,
                                                         const unsigned int& frameInterval,
                                                         const bool& prePersist,
                                                         const bool& postPersist,
                                                         GeomFrame** geomFrames,
                                                         spaceCore::Bvh<objCore::InterpolatedObj>** 
                                                            meshFrames){
    m_numberOfFrames = numberOfFrames;
    m_geomFrames = geomFrames;
    m_meshFrames = meshFrames;
    m_frameOffset = frameOffset;
    m_frameInterval = frameInterval;
    m_prePersist = prePersist;
    m_postPersist = postPersist;
}

HOST DEVICE AnimatedMeshContainer::AnimatedMeshContainer(AnimatedMeshContainerData data){
    m_numberOfFrames = data.m_numberOfFrames;
    m_geomFrames = data.m_geomFrames;
    m_meshFrames = data.m_meshFrames;
    m_id = data.m_id;
}

HOST DEVICE AnimatedMeshContainer::~AnimatedMeshContainer(){

}

HOST DEVICE GeomType AnimatedMeshContainer::GetType(){
    return MESH;
}

HOST DEVICE unsigned int AnimatedMeshContainer::GetID(){
    return m_id;
}

HOST DEVICE void AnimatedMeshContainer::Intersect(const rayCore::Ray& r, 
                                                  spaceCore::TraverseAccumulator& result){
    //translate frame into local frame space
    float clampedFrame = (r.m_frame - float(m_frameOffset)) / float(m_frameInterval);
    if((clampedFrame<0.0f && m_prePersist==false) || 
       (clampedFrame>=m_numberOfFrames && m_postPersist==false)){
        return;
    }
    clampedFrame = glm::clamp(clampedFrame, 0.0f, float(m_numberOfFrames-1));
    //ceil/floor to get frame indices
    unsigned int upperFrame = glm::ceil(clampedFrame);
    unsigned int lowerFrame = glm::floor(clampedFrame);
    float lerpWeight = clampedFrame - float(lowerFrame);
    //grab relevant transforms and LERP to get interpolated transform, apply inverse to ray
    glm::vec3 interpT = m_geomFrames[lowerFrame]->m_translation * (1.0f-lerpWeight) + 
                        m_geomFrames[upperFrame]->m_translation * lerpWeight;
    glm::vec3 interpR = m_geomFrames[lowerFrame]->m_rotation * (1.0f-lerpWeight) + 
                        m_geomFrames[upperFrame]->m_rotation * lerpWeight;
    glm::vec3 interpS = m_geomFrames[lowerFrame]->m_scale * (1.0f-lerpWeight) + 
                        m_geomFrames[upperFrame]->m_scale * lerpWeight;
    rayCore::Ray transformedR = r.Transform(utilityCore::buildInverseTransformationMatrix(interpT,
                                                                                          interpR,
                                                                                          interpS));
    //run intersection, transform result back into worldspace
    m_meshFrames[lowerFrame]->Traverse(transformedR, result);
    result.Transform(utilityCore::buildTransformationMatrix(interpT, interpR, interpS));
}
}
