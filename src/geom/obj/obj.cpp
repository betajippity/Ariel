// ObjCore4.0: An (improved) obj mesh wrangling library. Part of TAKUA Render.
// Written by Yining Karl Li
//
// File: obj.cpp
// Implements obj.hpp

#include <tbb/tbb.h>
#include <vector>
#include "obj.hpp"

namespace objCore {

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
    std::cout << m_numberOfVertices << " vertices" << std::endl;
    std::cout << m_numberOfNormals << " normals" << std::endl;
    std::cout << m_numberOfUVs << " uvs" << std::endl;
    std::cout << m_numberOfPolys << " polys" << std::endl;

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
                if(fv.w>=0){
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
                if(fv.w>=0){
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
}
