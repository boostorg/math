//  Copyright John Maddock 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include <iostream>
#include <iomanip>
#include <vector>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/relative_difference.hpp>
#include <boost/array.hpp>
#include "cuda_managed_ptr.hpp"
#include "stopwatch.hpp"

// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>

typedef float float_type;

/**
 * CUDA Kernel Device code
 *
 */
__global__ void cuda_test(const float_type *in, float_type *out, int numElements)
{
    using std::cos;
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        out[i] = boost::math::lgamma(in[i]);
    }
}

template <class T> struct table_type { typedef T type; };
typedef float_type T;
#define SC_(x) static_cast<T>(x)

#include "../test_gamma_data.ipp"

/**
 * Host main routine
 */
int main(void)
{
  try{
    // Consolidate the test data:
    std::vector<float_type> v;

    for(unsigned i = 0; i < factorials.size(); ++i)
       v.push_back(factorials[i][0]);
    for(unsigned i = 0; i < near_1.size(); ++i)
       v.push_back(near_1[i][0]);
    for(unsigned i = 0; i < near_2.size(); ++i)
       v.push_back(near_2[i][0]);
    for(unsigned i = 0; i < near_0.size(); ++i)
       v.push_back(near_0[i][0]);
    for(unsigned i = 0; i < near_m10.size(); ++i)
       v.push_back(near_m10[i][0]);
    for(unsigned i = 0; i < near_m55.size(); ++i)
       v.push_back(near_m55[i][0]);
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    // Print the vector length to be used, and compute its size
    int numElements = 50000;
    std::cout << "[Vector operation on " << numElements << " elements]" << std::endl;

    // Allocate the managed input vector A
    cuda_managed_ptr<float_type> input_vector(numElements);

    // Allocate the managed output vector C
    cuda_managed_ptr<float_type> output_vector(numElements);

    // Initialize the input vectors
    for (int i = 0; i < numElements; ++i)
    {
        int table_id = i % v.size();
        input_vector[i] = v[table_id];
    }

    // Launch the Vector Add CUDA Kernel
    int threadsPerBlock = 1024;
    int blocksPerGrid =(numElements + threadsPerBlock - 1) / threadsPerBlock;
    std::cout << "CUDA kernel launch with " << blocksPerGrid << " blocks of " << threadsPerBlock << " threads" << std::endl;
    
    watch w;
    cuda_test<<<blocksPerGrid, threadsPerBlock>>>(input_vector.get(), output_vector.get(), numElements);
    cudaDeviceSynchronize();
    std::cout << "CUDA kernal done in " << w.elapsed() << "s" << std::endl;
    
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        std::cerr << "Failed to launch vectorAdd kernel (error code " << cudaGetErrorString(err) << ")!" << std::endl;
        return EXIT_FAILURE;
    }

    // Verify that the result vector is correct
    std::vector<float_type> results;
    results.reserve(numElements);
    w.reset();
    for(int i = 0; i < numElements; ++i)
       results.push_back(boost::math::lgamma(input_vector[i]));
    double t = w.elapsed();
    // check the results
    for(int i = 0; i < numElements; ++i)
    {
        if (boost::math::epsilon_difference(output_vector[i], results[i]) > 300)
        {
            std::cerr << "Result verification failed at element " << i << "!" << std::endl;
            std::cerr << "Error rate was: " << boost::math::epsilon_difference(output_vector[i], results[i]) << "eps" << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cout << "Test PASSED with calculation time: " << t << "s" << std::endl;
    std::cout << "Done\n";
  }
  catch(const std::exception& e)
  {
    std::cerr << "Stopped with exception: " << e.what() << std::endl;
  }
  return 0;
}

