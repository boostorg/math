//  Copyright John Maddock 2016.
//  Copyright Matt Borland 2024 - 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Purposefully pull in headers that caused errors in the linked issue
#include <boost/math/special_functions.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/numeric_limits.hpp>
#include <boost/math/constants/constants.hpp>

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <random>
#include "cuda_managed_ptr.hpp"
#include "stopwatch.hpp"

// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>

typedef double float_type;

/**
 * CUDA Kernel Device code
 *
 */
__global__ void cuda_test(const float_type *in, float_type *out, int numElements)
{
    using std::cos;
    const int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        out[i] = cos(in[i]);
        if (out[i] > boost::math::numeric_limits<float_type>::max() || !boost::math::numeric_limits<float_type>::is_signed)
        {
            __trap();
        }
    }
}

/**
 * Host main routine
 */
int main()
{
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
    // Check some of our numeric_limits for viability
    std::mt19937_64 rng {42};
    std::uniform_real_distribution<float_type> dist(0, boost::math::constants::pi<float_type>());
    static_assert(boost::math::numeric_limits<float_type>::is_specialized, "Should be since it's a double");
    static_assert(boost::math::numeric_limits<float_type>::is_signed, "Should be since it's a double");

    for (int i = 0; i < numElements; ++i)
    {
        input_vector[i] = dist(rng);
    }

    // Launch the Vector Add CUDA Kernel
    int threadsPerBlock = 256;
    int blocksPerGrid =(numElements + threadsPerBlock - 1) / threadsPerBlock;
    std::cout << "CUDA kernel launch with " << blocksPerGrid << " blocks of " << threadsPerBlock << " threads" << std::endl;

    watch w;

    cuda_test<<<blocksPerGrid, threadsPerBlock>>>(input_vector.get(), output_vector.get(), numElements);
    cudaDeviceSynchronize();

    std::cout << "CUDA kernal done in: " << w.elapsed() << "s" << std::endl;

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
       results.push_back(std::cos(input_vector[i]));
    double t = w.elapsed();
    // check the results
    for(int i = 0; i < numElements; ++i)
    {
        if (boost::math::epsilon_difference(output_vector[i], results[i]) > 10)
        {
            std::cerr << "Result verification failed at element " << i << "!" << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cout << "Test PASSED, normal calculation time: " << t << "s" << std::endl;
    std::cout << "Done\n";

    return 0;
}
