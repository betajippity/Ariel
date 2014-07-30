// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: datastructures.hpp
// Some specialized common data structures that either STL doesn't have, or STL's version isn't
// quite what we need

#ifndef DATASTRUCTURES_HPP
#define DATASTRUCTURES_HPP

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#define SHARED __shared__
#else
#define HOST
#define DEVICE
#define SHARED
#endif

#include <tbb/tbb.h>
#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_vector.h>

//====================================
// Class Declarations
//====================================

//constant memory short stack implementation, NOT THREAD SAFE
template <typename T> class ShortStack {
    public:
        HOST DEVICE ShortStack();
        HOST DEVICE ~ShortStack();

        HOST DEVICE void Push(const T& item);
        HOST DEVICE T Pop();
        HOST DEVICE unsigned int Size();
        HOST DEVICE bool Empty();
        HOST DEVICE bool Full();
    
    private:
        T      m_stack[30];
        int    m_currentIndex;
};

//wrapper around tbb::concurrent_queue to add multipush/multipop. Is thread safe.
template <typename T> class MultiQueue {
    public:
        MultiQueue();
        ~MultiQueue();

        void Push(const T& item);
        void Push(const std::vector<T>& items);
        T Pop();
        std::vector<T> Pop(const unsigned int& count);
        unsigned int Size();
        bool Empty();

    private:
        tbb::concurrent_queue<T>   m_queue;
}; 

//wrapper around a vector that tracks the last element accessed and returns the next. 
//Is thread safe.
template <typename T> class LoopVector {
    public:
        LoopVector();
        ~LoopVector();

        void PushBack(const T& item);
        T GetElement();

    private:
        tbb::atomic<unsigned int>   m_currentIndex;
        tbb::concurrent_vector<T>   m_vector;
        tbb::atomic<unsigned int>   m_size;
};

#include "multiqueue.inl"
#include "shortstack.inl"
#include "loopvector.inl"

#endif
