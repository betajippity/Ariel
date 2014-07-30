// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: ray.hpp. Adapted from Takua Render.
// Classes for ray stuff

#ifndef RAY_HPP
#define RAY_HPP

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#define SHARED __shared__
#else
#define HOST
#define DEVICE
#define SHARED
#endif

#include "../utilities/utilities.h"

namespace rayCore {

//====================================
// Class Declarations
//====================================
class Ray {
    public:
        HOST DEVICE Ray();
        HOST DEVICE Ray(const glm::vec3& origin, const glm::vec3& direction, const float& frame,
                        const unsigned int& trackingID);
        HOST DEVICE Ray(const glm::vec3& origin, const glm::vec3& direction, const float& frame);
        HOST DEVICE ~Ray();

        HOST DEVICE void SetContents(const glm::vec3& origin, const glm::vec3& direction,
                                     const float& frame, const unsigned int& trackingID);

        HOST DEVICE glm::vec3 GetPointAlongRay(const float& distance) const;
        HOST DEVICE Ray Transform(const glm::mat4& m) const;

        glm::vec3       m_origin;
        glm::vec3       m_direction;
        unsigned int    m_trackingID;
        float           m_frame;
};

class Intersection {
    public:
        HOST DEVICE Intersection();
        HOST DEVICE Intersection(const bool& hit, const glm::vec3& point, const glm::vec3& normal,
                                 const glm::vec2& uv, const unsigned int& objectID, 
                                 const unsigned int& primID);
        HOST DEVICE ~Intersection();

        HOST DEVICE void SetContents(const bool& hit, const glm::vec3& point, 
                                     const glm::vec3& normal, const glm::vec2& uv, 
                                     const unsigned int& objectID, const unsigned int& primID);

        HOST DEVICE Intersection CompareClosestAgainst(const Intersection& hit, 
                                                       const glm::vec3& point);
        HOST DEVICE Intersection Transform(const glm::mat4& m);

        glm::vec3       m_point;
        glm::vec3       m_normal;
        glm::vec2       m_uv;
        unsigned int    m_objectID;
        unsigned int    m_primID;
        bool            m_hit;
};
}

#endif
