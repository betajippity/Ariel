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
    m_geomTransforms = NULL;
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
                                         GeomTransform** geomTransforms,
                                         spaceCore::Bvh<objCore::Obj>** meshFrames){
    m_numberOfFrames = numberOfFrames;
    m_geomTransforms = geomTransforms;
    m_meshFrames = meshFrames;
    m_frameOffset = frameOffset;
    m_frameInterval = frameInterval;
    m_prePersist = prePersist;
    m_postPersist = postPersist;
}

HOST DEVICE MeshContainer::MeshContainer(MeshContainerData data){
    m_numberOfFrames = data.m_numberOfFrames;
    m_geomTransforms = data.m_geomTransforms;
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

HOST DEVICE spaceCore::Bvh<objCore::Obj>* MeshContainer::GetMeshFrame(const float& frame){
    float clampedFrame = (frame - float(m_frameOffset)) / float(m_frameInterval);
    clampedFrame = glm::clamp(clampedFrame, 0.0f, float(m_numberOfFrames-1));
    unsigned int lowerFrame = glm::floor(clampedFrame);
    return m_meshFrames[lowerFrame];
}

HOST DEVICE spaceCore::Aabb MeshContainer::GetAabb(const float& frame){
    glm::mat4 transform;
    glm::mat4 inverse;
    if(GetTransforms(frame, transform, inverse)==true){
        spaceCore::Aabb box = GetMeshFrame(frame)->m_nodes[1].m_bounds;
        return box.Transform(transform);        
    }else{
        return spaceCore::Aabb();
    }
}       

HOST DEVICE bool MeshContainer::IsInFrame(const float& frame){
    float clampedFrame = (frame - float(m_frameOffset)) / float(m_frameInterval);
    if((clampedFrame<0.0f && m_prePersist==false) || 
       (clampedFrame>=m_numberOfFrames && m_postPersist==false)){
        return false;
    }
    return true;
}

HOST DEVICE bool MeshContainer::GetTransforms(const float& frame, glm::mat4& transform,
                                              glm::mat4& inversetransform){
    //translate frame into local frame space
    float clampedFrame = (frame - float(m_frameOffset)) / float(m_frameInterval);
    if((clampedFrame<0.0f && m_prePersist==false) || 
       (clampedFrame>=m_numberOfFrames && m_postPersist==false)){
        return false;
    }   
    clampedFrame = glm::clamp(clampedFrame, 0.0f, float(m_numberOfFrames-1));
    //ceil/floor to get frame indices
    unsigned int upperFrame = glm::ceil(clampedFrame);
    unsigned int lowerFrame = glm::floor(clampedFrame);
    float lerpWeight = clampedFrame - float(lowerFrame);
    //grab relevant transforms and LERP to get interpolated transform, apply inverse to ray
    glm::vec3 interpT = m_geomTransforms[lowerFrame]->m_translation * (1.0f-lerpWeight) + 
                        m_geomTransforms[upperFrame]->m_translation * lerpWeight;
    glm::vec3 interpR = m_geomTransforms[lowerFrame]->m_rotation * (1.0f-lerpWeight) + 
                        m_geomTransforms[upperFrame]->m_rotation * lerpWeight;
    glm::vec3 interpS = m_geomTransforms[lowerFrame]->m_scale * (1.0f-lerpWeight) + 
                        m_geomTransforms[upperFrame]->m_scale * lerpWeight;
    inversetransform = utilityCore::buildInverseTransformationMatrix(interpT, interpR, interpS);
    transform = utilityCore::buildTransformationMatrix(interpT, interpR, interpS);
    return true;
}

HOST DEVICE void MeshContainer::Intersect(const rayCore::Ray& r, 
                                          spaceCore::TraverseAccumulator& result){
    glm::mat4 transform;
    glm::mat4 inversetransform;
    if(GetTransforms(r.m_frame, transform, inversetransform)==false){
        return;
    }
    rayCore::Ray transformedR = r.Transform(inversetransform);
    //run intersection, transform result back into worldspace
    GetMeshFrame(r.m_frame)->Traverse(transformedR, result);
    result.Transform(transform);
}

HOST DEVICE bool MeshContainer::IsDynamic(){
    if(m_prePersist==true && m_postPersist==true && m_numberOfFrames==1){
        return false;
    }
    return true;
}

//====================================
// AnimatedMeshContainer Class
//====================================

HOST DEVICE AnimatedMeshContainer::AnimatedMeshContainer(){
    m_meshFrames = NULL;
    m_geomTransforms = NULL;
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
                                                         GeomTransform** geomTransforms,
                                                         spaceCore::Bvh<objCore::InterpolatedObj>** 
                                                            meshFrames){
    m_numberOfFrames = numberOfFrames;
    m_geomTransforms = geomTransforms;
    m_meshFrames = meshFrames;
    m_frameOffset = frameOffset;
    m_frameInterval = frameInterval;
    m_prePersist = prePersist;
    m_postPersist = postPersist;
}

HOST DEVICE AnimatedMeshContainer::AnimatedMeshContainer(AnimatedMeshContainerData data){
    m_numberOfFrames = data.m_numberOfFrames;
    m_geomTransforms = data.m_geomTransforms;
    m_meshFrames = data.m_meshFrames;
    m_id = data.m_id;
}

HOST DEVICE AnimatedMeshContainer::~AnimatedMeshContainer(){

}

HOST DEVICE spaceCore::Aabb AnimatedMeshContainer::GetAabb(const float& frame){
    glm::mat4 transform;
    glm::mat4 inverse;
    if(GetTransforms(frame, transform, inverse)==true){
        spaceCore::Aabb box = GetMeshFrame(frame)->m_nodes[1].m_bounds;
        return box.Transform(transform);        
    }else{
        return spaceCore::Aabb();
    }
}

HOST DEVICE GeomType AnimatedMeshContainer::GetType(){
    return ANIMMESH;
}

HOST DEVICE unsigned int AnimatedMeshContainer::GetID(){
    return m_id;
}

HOST DEVICE bool AnimatedMeshContainer::IsInFrame(const float& frame){
    float clampedFrame = (frame - float(m_frameOffset)) / float(m_frameInterval);
    if((clampedFrame<0.0f && m_prePersist==false) || 
       (clampedFrame>=m_numberOfFrames && m_postPersist==false)){
        return false;
    }
    return true;
}

HOST DEVICE bool AnimatedMeshContainer::GetTransforms(const float& frame, glm::mat4& transform,
                                                      glm::mat4& inversetransform){
    //translate frame into local frame space
    float clampedFrame = (frame - float(m_frameOffset)) / float(m_frameInterval);
    if((clampedFrame<0.0f && m_prePersist==false) || 
       (clampedFrame>=m_numberOfFrames && m_postPersist==false)){
        return false;
    }   
    clampedFrame = glm::clamp(clampedFrame, 0.0f, float(m_numberOfFrames-1));
    //ceil/floor to get frame indices
    unsigned int upperFrame = glm::ceil(clampedFrame);
    unsigned int lowerFrame = glm::floor(clampedFrame);
    float lerpWeight = clampedFrame - float(lowerFrame);
    //grab relevant transforms and LERP to get interpolated transform, apply inverse to ray
    glm::vec3 interpT = m_geomTransforms[lowerFrame]->m_translation * (1.0f-lerpWeight) + 
                        m_geomTransforms[upperFrame]->m_translation * lerpWeight;
    glm::vec3 interpR = m_geomTransforms[lowerFrame]->m_rotation * (1.0f-lerpWeight) + 
                        m_geomTransforms[upperFrame]->m_rotation * lerpWeight;
    glm::vec3 interpS = m_geomTransforms[lowerFrame]->m_scale * (1.0f-lerpWeight) + 
                        m_geomTransforms[upperFrame]->m_scale * lerpWeight;
    inversetransform = utilityCore::buildInverseTransformationMatrix(interpT, interpR, interpS);
    transform = utilityCore::buildTransformationMatrix(interpT, interpR, interpS);
    return true;
}

HOST DEVICE float AnimatedMeshContainer::GetInterpolationWeight(const float& frame){
    //translate frame into local frame space
    float clampedFrame = (frame - float(m_frameOffset)) / float(m_frameInterval);
    if((clampedFrame<0.0f && m_prePersist==false) || 
       (clampedFrame>=m_numberOfFrames && m_postPersist==false)){
        return false;
    }   
    clampedFrame = glm::clamp(clampedFrame, 0.0f, float(m_numberOfFrames-1));
    //ceil/floor to get frame indices
    unsigned int upperFrame = glm::ceil(clampedFrame);
    unsigned int lowerFrame = glm::floor(clampedFrame);
    float lerpWeight = clampedFrame - float(lowerFrame);
    return lerpWeight;
}

HOST DEVICE spaceCore::Bvh<objCore::InterpolatedObj>* AnimatedMeshContainer::GetMeshFrame(
                                                                            const float& frame){
    float clampedFrame = (frame - float(m_frameOffset)) / float(m_frameInterval);
    clampedFrame = glm::clamp(clampedFrame, 0.0f, float(m_numberOfFrames-1));
    unsigned int lowerFrame = glm::floor(clampedFrame);
    return m_meshFrames[lowerFrame];
}

HOST DEVICE void AnimatedMeshContainer::Intersect(const rayCore::Ray& r, 
                                                  spaceCore::TraverseAccumulator& result){
    glm::mat4 transform;
    glm::mat4 inversetransform;
    if(GetTransforms(r.m_frame, transform, inversetransform)==false){
        return;
    }
    rayCore::Ray transformedR = r.Transform(inversetransform);
    //run intersection, transform result back into worldspace
    float clampedFrame = (r.m_frame - float(m_frameOffset)) / float(m_frameInterval);
    clampedFrame = glm::clamp(clampedFrame, 0.0f, float(m_numberOfFrames-1));
    unsigned int lowerFrame = glm::floor(clampedFrame);
    m_meshFrames[lowerFrame]->Traverse(transformedR, result);
    result.Transform(transform);
}

HOST DEVICE bool AnimatedMeshContainer::IsDynamic(){
    return true;
}
}
