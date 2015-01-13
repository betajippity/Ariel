// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: bvh.inl. Adapted from Takua Render.
// Implements bvh.hpp

#ifndef BVH_INL
#define BVH_INL

#include "bvh.hpp"
#include "../utilities/datastructures.hpp"

namespace spaceCore {
   
template <typename T> Bvh<T>::Bvh(T basegeom){
    m_nodes = NULL;
    m_numberOfNodes = 0;
    m_referenceIndices = NULL;
    m_numberOfReferenceIndices = 0;
    m_basegeom = basegeom;
}

template <typename T> Bvh<T>::Bvh(){
    m_nodes = NULL;
    m_numberOfNodes = 0;
    m_referenceIndices = NULL;
    m_numberOfReferenceIndices = 0;
}

template <typename T> Bvh<T>::~Bvh(){
    /*if(m_nodes!=NULL){
        delete [] m_nodes;
    }
    if(m_referenceIndices!=NULL){
        delete [] m_referenceIndices;
    }*/
}

template <typename T> void Bvh<T>::Traverse(const rayCore::Ray& r, TraverseAccumulator& result){
    ShortStack<BvhNode*> stack;
    stack.Push(&m_nodes[1]);
    BvhNode* current = &m_nodes[1];
    //check for root node hit/miss
    float distanceToCurrent = current->FastIntersectionTest(r);
    if(distanceToCurrent<-0.5f){
        return;
    }else{
        float distanceToCurrent = 10000000000000.0f;
        while(stack.Empty()==false){
            //traverse and push back to stack until leaf node is reached
            bool empty = false;
            while(current->IsLeaf()==false && empty==false){
                BvhNode* left = &m_nodes[current->m_left];
                BvhNode* right = &m_nodes[current->m_right];
                //find nearest and farthest child
                float distanceToLeft = left->FastIntersectionTest(r);
                float distanceToRight = right->FastIntersectionTest(r);
                //if ray intersects both children, push current node onto stack and set current
                //to nearest child
                if(distanceToLeft>-0.5f && distanceToRight>-0.5f){
                    stack.Push(current);
                    if(distanceToLeft<distanceToRight){
                        current = left;
                    }else{
                        current = right;
                    }
                //else if intersect only the left or right child
                }else if(distanceToLeft>-0.5f && distanceToRight<-0.5f){
                    current = left;
                }else if(distanceToRight>-0.5f && distanceToLeft<-0.5f){
                    current = right;
                //else if intersect neither children, abort loop and move onto next node in stack
                }else{
                    empty = true;
                }
            }
            if(current->IsLeaf()){
                for(unsigned int i=0; i<current->m_numberOfReferences; i++){
                    unsigned int primID = m_referenceIndices[current->m_referenceOffset+i];
                    rayCore::Intersection rhit = m_basegeom.IntersectElement(primID, r); 
                    if(rhit.m_hit){
                        //make sure hit is actually in front of origin
                        glm::vec3 n = glm::normalize(rhit.m_point-r.m_origin);
                        float degree = glm::acos(glm::dot(n, r.m_direction));
                        if(degree<(PI/2.0f) || degree!=degree){
                            result.RecordIntersection(rhit, current->m_nodeid);
                        }
                    }
                }
            }
            //if stack is not empty, pop current node and pick farthest child
            if(stack.Empty()==false){
                current = stack.Pop();
                BvhNode* left = &m_nodes[current->m_left];
                BvhNode* right = &m_nodes[current->m_right];
                //find nearest and farthest child
                float distanceToLeft = left->FastIntersectionTest(r);
                float distanceToRight = right->FastIntersectionTest(r);
                if(distanceToLeft>=distanceToRight){
                    current = left;
                }else{
                    current = right;
                }
            }
        }
        return;
    }
}

template <typename T> void Bvh<T>::BuildBvh(const unsigned int& maxDepth){
    //assemble aabb list
    unsigned int numberOfAabbs = m_basegeom.GetNumberOfElements();
    spaceCore::Aabb* aabbs = new spaceCore::Aabb[numberOfAabbs];
    for(unsigned int i=0; i<numberOfAabbs; i++){
        aabbs[i] = m_basegeom.GetElementAabb(i);
    }
    //create vectors to store nodes during construction
    std::vector<BvhNode> tree;
    tree.reserve(std::pow(2.0f,(int)maxDepth));
    unsigned int nodeCount = 0;
    unsigned int layerCount = 0;
    std::vector< std::vector<unsigned int> > currentRefs; 
    currentRefs.reserve(std::pow(2.0f,(int)maxDepth));
    std::vector<BvhNode*> currentLayer; 
    std::vector<unsigned int> finalRefList;
    //create null node to place in index 0. tree begins at index 1.
    tree.push_back(BvhNode());
    nodeCount++;
//    BvhNode* nullnode = &tree[0];
    //create root node with aabb that is combined from all input aabbs
    tree.push_back(BvhNode());
    nodeCount++;
    BvhNode* root = &tree[1];
    for(unsigned int i=0; i<numberOfAabbs; i++){
        root->m_bounds.ExpandAabb(aabbs[i].m_min, aabbs[i].m_max);
    }
    //populate root node with all ids
    root->m_numberOfReferences = numberOfAabbs;
    root->m_referenceOffset = 0;
    root->m_nodeid = 1; 
    std::vector<unsigned int> rootreferences;
    rootreferences.reserve(numberOfAabbs);
    for(unsigned int i=0; i<numberOfAabbs; i++){
        rootreferences.push_back(aabbs[i].m_id);
    }
    currentLayer.push_back(root);
    currentRefs.push_back(rootreferences);
    //for each layer, split nodes and move on to next layer
    for(unsigned int layer=0; layer<maxDepth; layer++){
        //set up temp storage buffers for next layer
        std::vector< std::vector<unsigned int> > nextLayerRefs;
        std::vector<BvhNode*> nextLayer;
        //for each node in current layer, split and push result 
        unsigned int numberOfNodesInLayer = currentLayer.size(); 
        if(numberOfNodesInLayer>0){
            layerCount++;
        }
        for(unsigned int i=0; i<numberOfNodesInLayer; i++){
            //if <20 prims or last layer, make node a leaf node
            unsigned int refCount = currentRefs[i].size();
            BvhNode* node = currentLayer[i];
            if(refCount <= 5 || layer==maxDepth-1){
                node->m_left = 0;
                node->m_right = 0;
                node->m_numberOfReferences = refCount;
                node->m_referenceOffset = finalRefList.size();
                for(unsigned int j=0; j<node->m_numberOfReferences; j++){
                    finalRefList.push_back(currentRefs[i][j]);
                }   
            }else{ //else, find the best split via SAH and add children to next layer
                //split axis should be the longest axis
                Axis splitAxis = FindLongestAxis(node->m_bounds);
                //call SAH and get best split candidate
                unsigned int leftCount = 0;
                unsigned int rightCount = 0;
                float bestSplit = FindBestSplit(*node, currentRefs[i], splitAxis, 10, aabbs,
                                                leftCount, rightCount);
                if(leftCount==0 || rightCount==0){
                    leftCount = 0;
                    rightCount = 0;
                    if(splitAxis==axis_x){
                        splitAxis = axis_y;
                    }else if(splitAxis==axis_y){
                        splitAxis = axis_z;
                    }else if(splitAxis==axis_z){
                        splitAxis = axis_x;
                    }
                    bestSplit = FindBestSplit(*node, currentRefs[i], splitAxis, 10, aabbs,
                                                    leftCount, rightCount);
                }
                if(leftCount==0 || rightCount==0){
                    leftCount = 0;
                    rightCount = 0;
                    if(splitAxis==axis_x){
                        splitAxis = axis_y;
                    }else if(splitAxis==axis_y){
                        splitAxis = axis_z;
                    }else if(splitAxis==axis_z){
                        splitAxis = axis_x;
                    }
                    bestSplit = FindBestSplit(*node, currentRefs[i], splitAxis, 10, aabbs,
                                              leftCount, rightCount);
                }
                if(leftCount==0 || rightCount==0){
                    node->m_left = 0;
                    node->m_right = 0;
                    node->m_numberOfReferences = refCount;
                    node->m_referenceOffset = finalRefList.size();
                    for(unsigned int j=0; j<node->m_numberOfReferences; j++){
                        finalRefList.push_back(currentRefs[i][j]);
                    }   
                }else{
                    //create left and right nodes
                    tree.push_back(BvhNode());
                    nodeCount++;
                    unsigned int leftID = tree.size()-1;
                    BvhNode* left = &tree[leftID];
                    left->m_nodeid = leftID;
                    tree.push_back(BvhNode());
                    nodeCount++;
                    unsigned int rightID = tree.size()-1;
                    BvhNode* right = &tree[rightID];
                    right->m_nodeid = rightID;
                    node->m_left = leftID;
                    node->m_right = rightID;
                    //create lists of left and right objects based on centroid location
                    std::vector<unsigned int> leftRefs;
                    std::vector<unsigned int> rightRefs;
                    leftRefs.reserve(leftCount);
                    rightRefs.reserve(rightCount);
                    left->m_numberOfReferences = leftCount;
                    right->m_numberOfReferences = rightCount;
                    for(unsigned int j=0; j<node->m_numberOfReferences; j++){
                        unsigned int referenceIndex = currentRefs[i][j];
                        Aabb reference = aabbs[referenceIndex];
                        if(reference.m_centroid[splitAxis]<=bestSplit){
                            leftRefs.push_back(referenceIndex);
                            left->m_bounds.ExpandAabb(reference.m_min, reference.m_max);
                        }else{
                            rightRefs.push_back(referenceIndex);
                            right->m_bounds.ExpandAabb(reference.m_min, reference.m_max);
                        }
                    }
                    //add left and right children to next later
                    nextLayer.push_back(left);
                    nextLayerRefs.push_back(leftRefs);
                    nextLayer.push_back(right); 
                    nextLayerRefs.push_back(rightRefs);
                    //cleanup
                    leftRefs.clear();
                    rightRefs.clear();
                }
            }    
        }
        //clear current layer's temp buffers, swap in next layer
        currentLayer.clear();
        unsigned int numberOfCurrentRefs = currentRefs.size();
        for(unsigned int i=0; i<numberOfCurrentRefs; i++){
            currentRefs[i].clear();
        }
        currentRefs.clear();
        currentLayer = nextLayer;
        currentRefs = nextLayerRefs;
    }
    m_numberOfNodes = tree.size();
    m_nodes = new BvhNode[m_numberOfNodes];
    copy(tree.begin(), tree.end(), m_nodes);
    
    m_numberOfReferenceIndices = finalRefList.size();
    m_referenceIndices = new unsigned int[m_numberOfReferenceIndices];
    copy(finalRefList.begin(), finalRefList.end(), m_referenceIndices);

    delete [] aabbs;

    std::cout << "Built BVH with " << nodeCount << " nodes and depth " << layerCount << std::endl;
}

template <typename T> float Bvh<T>::FindBestSplit(BvhNode& node, 
                                                  std::vector<unsigned int>& references,
                                                  const Axis& direction, 
                                                  const unsigned int& quality, Aabb* aabbs,
                                                  unsigned int& leftCount, 
                                                  unsigned int& rightCount){
    int splitAxis = 0;
    if(direction==axis_y){splitAxis = 1;}
    else if(direction==axis_z){splitAxis = 2;}
    //create N even candidates, where N = quality
    float* candidates = new float[quality];
    float axisLength = node.m_bounds.m_max[splitAxis]-node.m_bounds.m_min[splitAxis];
    float evenSplitSize = (axisLength)/float(quality+1);
    for(unsigned int i=0; i<quality; i++){ //create even splits
        float candidate = node.m_bounds.m_min[splitAxis] + (i+1)*evenSplitSize;
        candidates[i] = candidate;
    }
    //Calculate all split costs and keep the lowest cost split
    unsigned int left = 0; unsigned int right = 0;
    unsigned int currentLeft = 0; unsigned int currentRight = 0;
    double bestCost = CalculateSplitCost(candidates[0], node, references, direction, aabbs, 
                                         left, right);
    float bestSplit = candidates[0];    
    for(unsigned int i=1; i<quality; i++){
        double cost = CalculateSplitCost(candidates[i], node, references, direction, aabbs, 
                                         currentLeft, currentRight);
        if(cost<bestCost){
            bestCost = cost;
            bestSplit = candidates[i];
            left = currentLeft;
            right = currentRight;
        }
    }
    leftCount = left;
    rightCount = right;
    delete [] candidates;
    return bestSplit;
}

template <typename T> double Bvh<T>::CalculateSplitCost(const float& split, const BvhNode& node, 
                                                        std::vector<unsigned int>& references, 
                                                        const Axis& direction, Aabb* aabbs, 
                                                        unsigned int& leftCount, 
                                                        unsigned int& rightCount){
    int splitAxis = 0;
    if(direction==axis_y){splitAxis = 1;}
    else if(direction==axis_z){splitAxis = 2;}
    Aabb leftBox; unsigned int left = 0;
    Aabb rightBox; unsigned int right = 0;
    //run through all references in current node and merge into either left or right
    for(unsigned int i=0; i<node.m_numberOfReferences; i++){
        unsigned int referenceIndex = references[i];
        Aabb reference = aabbs[referenceIndex];
        if(reference.m_centroid[splitAxis]<=split){
            left++;
            leftBox.ExpandAabb(reference.m_min, reference.m_max);
        }else{
            right++;
            rightBox.ExpandAabb(reference.m_min, reference.m_max);
        }
    }
    double leftCost = left * leftBox.CalculateSurfaceArea();
    double rightCost = right * rightBox.CalculateSurfaceArea();
    leftCount = left; rightCount = right;
    // std::cout << left << " "  << right << " " << fabs(leftCost-rightCost) << std::endl;
    if(left==0 || right==0){
        return pow(10.0f,100);
    }else{
        return fabs(leftCost-rightCost);
    }
}

template <typename T> Axis Bvh<T>::FindLongestAxis(const Aabb& box){
    //find max axis and return. pretty simple.
    float axisXlength = box.m_max.x - box.m_min.x;
    float axisYlength = box.m_max.y - box.m_min.y; 
    float axisZlength = box.m_max.z - box.m_min.z;
    if(axisXlength>=axisYlength && axisXlength>=axisZlength){
        return axis_x;
    }else if(axisYlength>=axisXlength && axisYlength>=axisZlength){
        return axis_y;
    }else{
        return axis_z;
    }   
}
}

#endif
