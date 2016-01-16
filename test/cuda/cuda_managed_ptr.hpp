//  Copyright John Maddock 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUDA_MANAGED_PTR_HPP
#define BOOST_MATH_CUDA_MANAGED_PTR_HPP

#ifdef _MSC_VER
#pragma once
#endif

class managed_holder_base
{
protected:
   static int count;
   managed_holder_base() { ++count; }
   ~managed_holder_base()
   {
      if(0 == --count)
         cudaDeviceSynchronize();
   }
};

int managed_holder_base::count = 0;

template <class T>
class cuda_managed_ptr
{
   T* data;
   cuda_managed_ptr(const cuda_managed_ptr&) = delete;
   cuda_managed_ptr& operator=(cuda_managed_ptr const&) = delete;
   void free()
   {
      if(data)
      {
         cudaDeviceSynchronize();
         cudaFree(data);
      }
   }
public:
   cuda_managed_ptr() : data(0) {}
   cuda_managed_ptr(std::size_t n)
   {
      cudaError_t err = cudaSuccess;
      void *ptr;
      err = cudaMallocManaged(&ptr, n * sizeof(T));
      if(err != cudaSuccess)
         throw std::runtime_error(cudaGetErrorString(err));
      cudaDeviceSynchronize();
      data = static_cast<T*>(ptr);
   }
   cuda_managed_ptr(cuda_managed_ptr&& o)
   {
      data = o.data;
      o.data = 0;
   }
   cuda_managed_ptr& operator=(cuda_managed_ptr&& o)
   {
      free();
      data = o.data;
      o.data = 0;
      return *this;
   }
   ~cuda_managed_ptr()
   {
      free();
   }

   class managed_holder : managed_holder_base
   {
      T* pdata;
   public:
      managed_holder(T* p) : managed_holder_base(), pdata(p) {}
      managed_holder(const managed_holder& o) : managed_holder_base(), pdata(o.pdata) {}
      operator T* () { return pdata; }
      T& operator[] (std::size_t n) { return pdata[n]; }
   };
   class const_managed_holder : managed_holder_base
   {
      const T* pdata;
   public:
      const_managed_holder(T* p) : managed_holder_base(), pdata(p) {}
      const_managed_holder(const managed_holder& o) : managed_holder_base(), pdata(o.pdata) {}
      operator const T* () { return pdata; }
      const T& operator[] (std::size_t n) { return pdata[n]; }
   };

   managed_holder get() { return managed_holder(data); }
   const_managed_holder get()const { return data; }
   T& operator[](std::size_t n) { return data[n]; }
   const T& operator[](std::size_t n)const { return data[n]; }
};

#endif


