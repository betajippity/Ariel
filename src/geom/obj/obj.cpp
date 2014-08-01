// ObjCore4.0: An (improved) obj mesh wrangling library. Part of TAKUA Render.
// Written by Yining Karl Li
//
// File: obj.cpp
// Implements obj.hpp

#include <tbb/tbb.h>
#include <vector>
#include "obj.hpp"

namespace objCore {

//====================================
// Obj Class
//====================================

Obj::Obj(){
    m_numberOfVertices = 0;
    m_numberOfNormals = 0;
    m_numberOfUVs = 0;
    m_numberOfPolys = 0;
    m_vertices = NULL;
    m_normals = NULL;
    m_uvs = NULL;
    m_polyVertexIndices = NULL;
    m_polyNormalIndices = NULL;
    m_polyUVIndices = NULL;
    m_keep = false;
}

Obj::Obj(const std::string& filename){
    ReadObj(filename);
    m_keep = false;
}

Obj::~Obj(){
    if(!m_keep){ //if keep flag is set true, contents won't be destroyed. Used for memory transfer.
        delete [] m_vertices;
        delete [] m_polyVertexIndices;
        if(m_numberOfNormals>0){
            delete [] m_normals;
            delete [] m_polyNormalIndices;
        }
        if(m_numberOfUVs){
            delete [] m_uvs;
            delete [] m_polyUVIndices;
        }
    }
}

void Obj::BakeTransform(const glm::mat4& transform){
    tbb::parallel_for(tbb::blocked_range<unsigned int>(0,m_numberOfVertices),
        [=](const tbb::blocked_range<unsigned int>& r){
            for(unsigned int i=r.begin(); i!=r.end(); ++i){ 
                m_vertices[i] = glm::vec3(transform * glm::vec4(m_vertices[i], 1.0f));
            }
        }
    );
    tbb::parallel_for(tbb::blocked_range<unsigned int>(0,m_numberOfNormals),
        [=](const tbb::blocked_range<unsigned int>& r){
            for(unsigned int i=r.begin(); i!=r.end(); ++i){ 
                m_normals[i] = glm::normalize(glm::vec3(transform * glm::vec4(m_normals[i], 0.0f)));
            }
        }
    );
}

bool Obj::ReadObj(const std::string& filename){
    PrereadObj(filename);

    m_vertices = new glm::vec3[m_numberOfVertices];
    m_normals = new glm::vec3[m_numberOfNormals];
    m_uvs = new glm::vec2[m_numberOfUVs];
    m_polyVertexIndices = new glm::uvec4[m_numberOfPolys];
    m_polyNormalIndices = new glm::uvec4[m_numberOfPolys];
    m_polyUVIndices = new glm::uvec4[m_numberOfPolys];

    unsigned int currentVertex = 0;
    unsigned int currentNormal = 0;
    unsigned int currentUV = 0;
    unsigned int currentFace = 0;

    //Load loop
    std::ifstream fp_in;
    char* fname = (char*)filename.c_str();
    fp_in.open(fname);
    if(fp_in.is_open()){
        while(fp_in.good()){
            std::string line;
            getline(fp_in, line);
            if(!line.empty()){
                std::vector<std::string> tokens = utilityCore::tokenizeString(line, " ");
                if(tokens.size()>1){
                    if(tokens[0].compare("v")==0){
                        m_vertices[currentVertex] = glm::vec3(atof(tokens[1].c_str()),
                                                              atof(tokens[2].c_str()),
                                                              atof(tokens[3].c_str()));
                        currentVertex++;
                    }else if(tokens[0].compare("vn")==0){
                        m_normals[currentNormal] = glm::vec3(atof(tokens[1].c_str()),
                                                             atof(tokens[2].c_str()),
                                                             atof(tokens[3].c_str()));
                        currentNormal++;
                    }else if(tokens[0].compare("vt")==0){
                        m_uvs[currentUV] = glm::vec2(atof(tokens[1].c_str()),
                                                     atof(tokens[2].c_str()));
                        currentUV++;
                    }else if(tokens[0].compare("f")==0){
                        glm::vec4 vertices(0);
                        glm::vec4 normals(0);
                        glm::vec4 uvs(0);
                        unsigned int loops = std::min(4, (int)tokens.size()-1)+1;
                        for(unsigned int i=1; i<loops; i++){
                            std::vector<std::string> faceindices =
                                                     utilityCore::tokenizeString(tokens[i], "/");
                            vertices[i-1] = atoi(faceindices[0].c_str());
                            if(faceindices.size()==2 && m_numberOfNormals>1){
                                normals[i-1] = atoi(faceindices[1].c_str());
                            }else if(faceindices.size()==2 && m_numberOfUVs>1){
                                uvs[i-1] = atoi(faceindices[1].c_str());
                            }else if(faceindices.size()==3){
                                uvs[i-1] = atoi(faceindices[1].c_str());
                                normals[i-1] = atoi(faceindices[2].c_str());
                            }
                        }
                        m_polyVertexIndices[currentFace] = vertices;
                        m_polyNormalIndices[currentFace] = normals;
                        m_polyUVIndices[currentFace] = uvs;
                        currentFace++;
                    }
                }               
            }
        }
    } 
    fp_in.close();

    if(m_numberOfUVs==0){
        std::cout << "No UVs found, creating default UVs..." << std::endl;
        delete [] m_uvs;
        m_uvs = new glm::vec2[4];
        m_uvs[0] = glm::vec2(0,0);
        m_uvs[1] = glm::vec2(0,1);
        m_uvs[2] = glm::vec2(1,1);
        m_uvs[3] = glm::vec2(1,0);
        m_numberOfUVs = 4;
        for(unsigned int i=0; i<m_numberOfPolys; i++){
            m_polyUVIndices[i] = glm::uvec4(1,2,3,4);
        }
    }

    if(m_numberOfNormals==0){
        std::cout << "No normals found, creating default normals..." << std::endl;
        delete [] m_normals;
        m_numberOfNormals = m_numberOfPolys;
        m_normals = new glm::vec3[m_numberOfPolys];
        for(unsigned int i=0; i<m_numberOfPolys; i++){
            glm::uvec4 vertexIndices = m_polyVertexIndices[i];
            glm::vec3 v0 = m_vertices[vertexIndices[0]-1];
            glm::vec3 v1 = m_vertices[vertexIndices[1]-1];
            glm::vec3 v2 = m_vertices[vertexIndices[2]-1];
            glm::vec3 n0 = v0-v1;
            glm::vec3 n1 = v2-v1;
            m_normals[i] = glm::normalize(glm::cross(n1,n0));
            m_polyNormalIndices[i] = glm::uvec4(i+1,i+1,i+1,i+1);
        }
    }

    std::cout << "Read obj from " << filename << std::endl;
    // std::cout << m_numberOfVertices << " vertices" << std::endl;
    // std::cout << m_numberOfNormals << " normals" << std::endl;
    // std::cout << m_numberOfUVs << " uvs" << std::endl;
    // std::cout << m_numberOfPolys << " polys" << std::endl;

    return true;
}

bool Obj::WriteObj(const std::string& filename){
    std::ofstream outputFile(filename.c_str());
    if(outputFile.is_open()){
        //write out vertices
        for(unsigned int i=0; i<m_numberOfVertices; i++){
            glm::vec3 v = m_vertices[i];
            outputFile << "v " << v.x << " " << v.y << " " << v.z << "\n";
        }
        //write out uvs
        for(unsigned int i=0; i<m_numberOfUVs; i++){
            glm::vec2 uv = m_uvs[i];
            outputFile << "vt " << uv.x << " " << uv.y << "\n";
        }
        //write out normals
        for(unsigned int i=0; i<m_numberOfNormals; i++){
            glm::vec3 n = m_normals[i];
            outputFile << "vn " << n.x << " " << n.y << " " << n.z << "\n";
        }
        //Write out faces
        if(m_numberOfUVs==0 && m_numberOfNormals==0){
            for(unsigned int i=0; i<m_numberOfPolys; i++){
                glm::uvec4 fv = m_polyVertexIndices[i];
                outputFile << "f ";
                outputFile << fv.x << " " << fv.y << " " << fv.z << " ";
                if(fv.w!=0){
                    outputFile << fv.w;
                }
                outputFile << "\n";
            }   
        }else{
            for(unsigned int i=0; i<m_numberOfPolys; i++){
                glm::uvec4 fv = m_polyVertexIndices[i];
                glm::uvec4 fn = m_polyNormalIndices[i];
                glm::uvec4 fuv = m_polyUVIndices[i];
                outputFile << "f ";
                outputFile << fv.x << "/" << fuv.x << "/" << fn.x << " ";
                outputFile << fv.y << "/" << fuv.y << "/" << fn.y << " ";
                outputFile << fv.z << "/" << fuv.z << "/" << fn.z << " ";
                if(fv.w!=0){
                    outputFile << fv.w << "/" << fuv.w << "/" << fn.w << " ";
                }
                outputFile << "\n";
            }            
        }
        outputFile.close();
        std::cout << "Wrote obj file to " << filename << std::endl;
        return true;
    }else{
        std::cout << "Error: Unable to write to " << filename << std::endl;
        return false;
    }
}

void Obj::PrereadObj(const std::string& filename){
    std::ifstream fp_in;
    unsigned int vertexCount = 0;
    unsigned int normalCount = 0;
    unsigned int uvCount = 0;
    unsigned int faceCount = 0;
    char* fname = (char*)filename.c_str();
    fp_in.open(fname);
    if(fp_in.is_open()){
        while(fp_in.good()){
            std::string line;
            getline(fp_in, line);
            if(!line.empty()){
                std::string header = utilityCore::getFirstNCharactersOfString(line, 2);
                if(header.compare("v ")==0){
                    vertexCount++;
                }else if(header.compare("vt")==0){
                    uvCount++;
                }else if(header.compare("vn")==0){
                    normalCount++;
                }else if(header.compare("f ")==0){
                    faceCount++;
                }
            }
        }
    } 
    fp_in.close();
    m_numberOfVertices = vertexCount;
    m_numberOfNormals = normalCount;
    m_numberOfUVs = uvCount;
    m_numberOfPolys = faceCount;
}

/*Return the requested face from the mesh, unless the index is out of range, 
in which case return a face of area zero*/
HOST DEVICE Poly Obj::GetPoly(const unsigned int& polyIndex){
    Point pNull(glm::vec3(0,0,0), glm::vec3(0,1,0), glm::vec2(0,0));
    if(polyIndex>=m_numberOfPolys){
        return Poly(pNull, pNull, pNull, polyIndex);
    }else{
        glm::uvec4 vertexIndices = m_polyVertexIndices[polyIndex];
        glm::uvec4 normalIndices = m_polyNormalIndices[polyIndex];
        glm::uvec4 uvIndices = m_polyUVIndices[polyIndex];
        Point p1(m_vertices[vertexIndices.x-1], m_normals[normalIndices.x-1], m_uvs[uvIndices.x-1]);
        Point p2(m_vertices[vertexIndices.y-1], m_normals[normalIndices.y-1], m_uvs[uvIndices.y-1]);
        Point p3(m_vertices[vertexIndices.z-1], m_normals[normalIndices.z-1], m_uvs[uvIndices.z-1]);
        if(vertexIndices.w==0){
            return Poly(p1, p2, p3, p1, polyIndex);
        }else{
            Point p4(m_vertices[vertexIndices.w-1], m_normals[normalIndices.w-1], 
                     m_uvs[uvIndices.w-1]);
            return Poly(p1, p2, p3, p4, polyIndex);
        }
        
    }
}

//Applies given transforglm::mation glm::matrix to the given point
HOST DEVICE Point Obj::TransformPoint(const Point& p, const glm::mat4& m){
    Point r = p;
    r.m_position = glm::vec3( utilityCore::multiply(m,glm::vec4(p.m_position,1.0)) );
    r.m_normal = glm::normalize(glm::vec3( utilityCore::multiply(m,glm::vec4(p.m_normal,0.0)) ));
    return r;
}

//Applies given transforglm::mation glm::matrix to the given poly
HOST DEVICE Poly Obj::TransformPoly(const Poly& p, const glm::mat4& m){
    Poly r;
    r.m_vertex0 = TransformPoint(p.m_vertex0, m);
    r.m_vertex1 = TransformPoint(p.m_vertex1, m);
    r.m_vertex2 = TransformPoint(p.m_vertex2, m);
    r.m_vertex3 = TransformPoint(p.m_vertex3, m);
    r.m_id = p.m_id;
    return r;
}

HOST DEVICE unsigned int Obj::GetNumberOfElements(){
    return m_numberOfPolys;
}

HOST DEVICE spaceCore::Aabb Obj::GetElementAabb(const unsigned int& primID){
    glm::uvec4 vertexIndices = m_polyVertexIndices[primID];
    glm::vec3 v0 = m_vertices[vertexIndices.x-1];
    glm::vec3 v1 = m_vertices[vertexIndices.y-1];
    glm::vec3 v2 = m_vertices[vertexIndices.z-1];
    glm::vec3 v3 = v0;
    if(vertexIndices.w>0){
        v3 = m_vertices[vertexIndices.w-1];  
    }
    glm::vec3 min = glm::min(glm::min(glm::min(v0, v1), v2), v3);
    glm::vec3 max = glm::max(glm::max(glm::max(v0, v1), v2), v3);
    //if v0 and v3 are the same, it's a triangle! else, it's a quad. TODO: find better way to
    //handle this check
    glm::vec3 centroid = glm::vec3(0.0f);
    if(vertexIndices.w>0){
        centroid = (v0 + v1 + v2)/3.0f;
    }else{
        centroid = (v0 + v1 + v2 + v3)/4.0f;
    }
    return spaceCore::Aabb(min, max, centroid, primID);
}

HOST DEVICE rayCore::Intersection Obj::IntersectElement(const unsigned int& primID, 
                                                        const rayCore::Ray& r){
    //check if quad by comparing x and w components
    glm::uvec4 vi = m_polyVertexIndices[primID];  
    //triangle case
    if(vi.w==0){
        return TriangleTest(primID, r, false);
    }else{ //quad case
        return QuadTest(primID, r); 
    }
}

HOST DEVICE rayCore::Intersection Obj::QuadTest(const unsigned int& polyIndex, 
                                                const rayCore::Ray& r){
    rayCore::Intersection intersect = TriangleTest(polyIndex, r, false);
    if(intersect.m_hit){
         return intersect;
    }else{
        return TriangleTest(polyIndex, r, true);
    }
}

HOST DEVICE inline rayCore::Intersection Obj::RayTriangleTest(const glm::vec3& v0,
                                                              const glm::vec3& v1,
                                                              const glm::vec3& v2,
                                                              const glm::vec3& n0,
                                                              const glm::vec3& n1,
                                                              const glm::vec3& n2,
                                                              const glm::vec2& u0,
                                                              const glm::vec2& u1,
                                                              const glm::vec2& u2,
                                                              const rayCore::Ray& r){
    //grab points and edges from poly
    glm::vec3 edge1 = v1-v0;
    glm::vec3 edge2 = v2-v0;
    //calculate determinant
    glm::vec3 rdirection = r.m_direction;
    glm::vec3 pvec = glm::cross(rdirection, edge2);
    float det = glm::dot(edge1, pvec);
    if(det == 0.0f){
        return rayCore::Intersection();
    }else{
        float inv_det = 1.0f/det;
        glm::vec3 tvec = r.m_origin - v1;
        //calculate barycentric coord
        glm::vec3 bary;
        bary.x = glm::dot(tvec, pvec) * inv_det;
        glm::vec3 qvec = glm::cross(tvec, edge1);
        bary.y = glm::dot(rdirection, qvec) * inv_det;
        //calcualte distance from ray origin to intersection
        float t = glm::dot(edge2, qvec) * inv_det;
        bool hit = (bary.x >= 0.0f && bary.y >= 0.0f && (bary.x + bary.y) <= 1.0f);
        if(hit){
            bary.z = 1.0f - bary.x - bary.y;
            // glm::vec3 hitPoint = r.m_origin + t*r.m_direction;
            glm::vec3 hitPoint = r.GetPointAlongRay(t);
            glm::vec3 hitNormal = (n0 * bary.z)+
                                  (n1 * bary.x)+
                                  (n2 * bary.y);
            hitNormal = hitNormal/glm::length(hitNormal);
            glm::vec2 hitUV = (u0 * bary.z)+ 
                              (u1 * bary.x)+
                              (u2 * bary.y);
            return rayCore::Intersection(true, hitPoint, hitNormal, hitUV, 0, 0);
        }else{
            return rayCore::Intersection();
        }
    }
}

HOST DEVICE rayCore::Intersection Obj::TriangleTest(const unsigned int& polyIndex,
                                                    const rayCore::Ray& r,
                                                    const bool& checkQuad){
    //grab indices
    glm::uvec4 vi = m_polyVertexIndices[polyIndex];
    glm::uvec4 ni = m_polyNormalIndices[polyIndex];
    glm::uvec4 ui = m_polyUVIndices[polyIndex];
    //grab points, edges, uvs from poly
    glm::vec3 v0 = m_vertices[vi.x-1];
    glm::vec3 v1 = m_vertices[vi.y-1];
    glm::vec3 v2 = m_vertices[vi.z-1];
    glm::vec3 n0 = m_normals[ni.x-1];
    glm::vec3 n1 = m_normals[ni.y-1];
    glm::vec3 n2 = m_normals[ni.z-1];
    glm::vec2 u0 = m_uvs[ui.x-1];
    glm::vec2 u1 = m_uvs[ui.y-1];
    glm::vec2 u2 = m_uvs[ui.z-1];
    if(checkQuad==true){
        v1 = m_vertices[vi.w-1];
        n1 = m_normals[ni.w-1];
        u1 = m_uvs[ui.w-1];
    }
    rayCore::Intersection result = Obj::RayTriangleTest(v0, v1, v2, n0, n1, n2, u0, u1, u2, r);
    result.m_objectID = m_id;
    result.m_primID = polyIndex;
    return result;
}

//====================================
// InterpolatedObj Class
//====================================

InterpolatedObj::InterpolatedObj(){
    m_obj0 = NULL;
    m_obj1 = NULL;
}

/*Right now only prints a warning if objs have mismatched topology, must make this do something
better in the future since mismatched topology leads to Very Bad Things*/
InterpolatedObj::InterpolatedObj(objCore::Obj* obj0, objCore::Obj* obj1){
    m_obj0 = obj0;
    m_obj1 = obj1;
    if(obj0->m_numberOfPolys!=obj1->m_numberOfPolys){
        std::cout << "Warning: Attempted to create InterpolatedObj with mismatched topology!"
                  << std::endl;
    }
}

InterpolatedObj::~InterpolatedObj(){
}

HOST DEVICE rayCore::Intersection InterpolatedObj::IntersectElement(const unsigned int& primID, 
                                                                    const rayCore::Ray& r){
    //check if quad by comparing x and w components
    glm::uvec4 vi = m_obj0->m_polyVertexIndices[primID];  
    //triangle case
    if(vi.w==0){
        return TriangleTest(primID, r, false);
    }else{ //quad case
        return QuadTest(primID, r); 
    }
}

HOST DEVICE rayCore::Intersection InterpolatedObj::TriangleTest(const unsigned int& polyIndex, 
                                                                const rayCore::Ray& r,
                                                                const bool& checkQuad){
    //make sure interp value is between 0 and 1
    float clampedInterp = r.m_frame - glm::floor(r.m_frame); 
    //grab indices
    glm::uvec4 vi0 = m_obj0->m_polyVertexIndices[polyIndex];
    glm::uvec4 ni0 = m_obj0->m_polyNormalIndices[polyIndex];
    glm::uvec4 ui0 = m_obj0->m_polyUVIndices[polyIndex];
    glm::uvec4 vi1 = m_obj1->m_polyVertexIndices[polyIndex];
    glm::uvec4 ni1 = m_obj1->m_polyNormalIndices[polyIndex];
    glm::uvec4 ui1 = m_obj1->m_polyUVIndices[polyIndex];
    //grab points, edges, uvs from poly
    glm::vec3 v0 = m_obj0->m_vertices[vi0.x-1] * (1.0f-clampedInterp) + 
                   m_obj1->m_vertices[vi1.x-1] * clampedInterp;
    glm::vec3 v1 = m_obj0->m_vertices[vi0.y-1] * (1.0f-clampedInterp) + 
                   m_obj1->m_vertices[vi1.y-1] * clampedInterp;
    glm::vec3 v2 = m_obj0->m_vertices[vi0.z-1] * (1.0f-clampedInterp) + 
                   m_obj1->m_vertices[vi1.z-1] * clampedInterp;
    glm::vec3 n0 = m_obj0->m_normals[ni0.x-1] * (1.0f-clampedInterp) + 
                   m_obj1->m_normals[ni1.x-1] * clampedInterp;
    glm::vec3 n1 = m_obj0->m_normals[ni0.y-1] * (1.0f-clampedInterp) + 
                   m_obj1->m_normals[ni1.y-1] * clampedInterp;
    glm::vec3 n2 = m_obj0->m_normals[ni0.z-1] * (1.0f-clampedInterp) + 
                   m_obj1->m_normals[ni1.z-1] * clampedInterp;
    glm::vec2 u0 = m_obj0->m_uvs[ui0.x-1] * (1.0f-clampedInterp) + 
                   m_obj1->m_uvs[ui1.x-1] * clampedInterp;
    glm::vec2 u1 = m_obj0->m_uvs[ui0.y-1] * (1.0f-clampedInterp) + 
                   m_obj1->m_uvs[ui1.y-1] * clampedInterp;
    glm::vec2 u2 = m_obj0->m_uvs[ui0.z-1] * (1.0f-clampedInterp) + 
                   m_obj1->m_uvs[ui1.z-1] * clampedInterp;
    if(checkQuad==true){
        v1 = m_obj0->m_vertices[vi0.w-1] * (1.0f-clampedInterp) + 
             m_obj1->m_vertices[vi1.w-1] * clampedInterp;
        n1 = m_obj0->m_normals[ni0.w-1] * (1.0f-clampedInterp) + 
             m_obj1->m_normals[ni1.w-1] * clampedInterp;
        u1 = m_obj0->m_uvs[ui0.w-1] * (1.0f-clampedInterp) + 
             m_obj1->m_uvs[ui1.w-1] * clampedInterp;
    }
    rayCore::Intersection result = Obj::RayTriangleTest(v0, v1, v2, n0/glm::length(n0), 
                                                        n1/glm::length(n1), n2/glm::length(n2), 
                                                        u0, u1, u2, r);
    result.m_objectID = m_obj0->m_id;
    result.m_primID = polyIndex;
    return result;
}

HOST DEVICE rayCore::Intersection InterpolatedObj::QuadTest(const unsigned int& polyIndex,
                                                            const rayCore::Ray& r){
    rayCore::Intersection intersect = TriangleTest(polyIndex, r, false);
    if(intersect.m_hit){
         return intersect;
    }else{
        return TriangleTest(polyIndex, r, true);
    }
}

/*Calls GetPoly for both referenced objs and returns a single interpolated poly. */
HOST DEVICE Poly InterpolatedObj::GetPoly(const unsigned int& polyIndex, 
                                          const float& interpolation){
    Poly p0 = m_obj0->GetPoly(polyIndex);
    Poly p1 = m_obj1->GetPoly(polyIndex);
    //make sure interp value is between 0 and 1
    float clampedInterp = interpolation - glm::floor(interpolation); 
    Poly p;
    p.m_vertex0.m_position = p0.m_vertex0.m_position * (1.0f-clampedInterp) + 
                             p1.m_vertex0.m_position * clampedInterp;
    p.m_vertex1.m_position = p0.m_vertex1.m_position * (1.0f-clampedInterp) + 
                             p1.m_vertex1.m_position * clampedInterp;
    p.m_vertex2.m_position = p0.m_vertex2.m_position * (1.0f-clampedInterp) + 
                             p1.m_vertex2.m_position * clampedInterp;
    p.m_vertex3.m_position = p0.m_vertex3.m_position * (1.0f-clampedInterp) + 
                             p1.m_vertex3.m_position * clampedInterp;
    p.m_vertex0.m_normal = p0.m_vertex0.m_normal * (1.0f-clampedInterp) + 
                           p1.m_vertex0.m_normal * clampedInterp;
    p.m_vertex1.m_normal = p0.m_vertex1.m_normal * (1.0f-clampedInterp) + 
                           p1.m_vertex1.m_normal * clampedInterp;
    p.m_vertex2.m_normal = p0.m_vertex2.m_normal * (1.0f-clampedInterp) + 
                           p1.m_vertex2.m_normal * clampedInterp;
    p.m_vertex3.m_normal = p0.m_vertex3.m_normal * (1.0f-clampedInterp) + 
                           p1.m_vertex3.m_normal * clampedInterp;
    p.m_vertex0.m_uv = p0.m_vertex0.m_uv * (1.0f-clampedInterp) + 
                       p1.m_vertex0.m_uv * clampedInterp;
    p.m_vertex1.m_uv = p0.m_vertex1.m_uv * (1.0f-clampedInterp) + 
                       p1.m_vertex1.m_uv * clampedInterp;
    p.m_vertex2.m_uv = p0.m_vertex2.m_uv * (1.0f-clampedInterp) + 
                       p1.m_vertex2.m_uv * clampedInterp;
    p.m_vertex3.m_uv = p0.m_vertex3.m_uv * (1.0f-clampedInterp) + 
                       p1.m_vertex3.m_uv * clampedInterp;
    return p;
}

HOST DEVICE spaceCore::Aabb InterpolatedObj::GetElementAabb(const unsigned int& primID){
    spaceCore::Aabb aabb0 = m_obj0->GetElementAabb(primID);
    spaceCore::Aabb aabb1 = m_obj1->GetElementAabb(primID);
    glm::vec3 combinedMin = glm::min(aabb0.m_min, aabb1.m_min);
    glm::vec3 combinedMax = glm::max(aabb0.m_max, aabb1.m_max);
    glm::vec3 combinedCentroid = (aabb0.m_centroid + aabb1.m_centroid) / 2.0f;
    return spaceCore::Aabb(combinedMin, combinedMax, combinedCentroid, aabb0.m_id);
}

HOST DEVICE unsigned int InterpolatedObj::GetNumberOfElements(){
    return m_obj0->GetNumberOfElements();
}
}
