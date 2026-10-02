//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <random>
#include <vector>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/relative_difference.hpp>
#include "cuda_managed_ptr.hpp"
#include "stopwatch.hpp"

// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>

typedef float float_type;

/**
 * CUDA Kernel Device code
 *
 */
__global__ void cuda_test(const float_type *in1, const float_type *in2, float_type *out, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        out[i] = boost::math::ligamma(in1[i], in2[i]);
    }
}

/**
 * Host main routine
 */
int main(void)
{
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    // Print the vector length to be used, and compute its size
    int numElements = 50000;
    std::cout << "[Vector operation on " << numElements << " elements]" << std::endl;

    // Allocate the managed input vectors
    cuda_managed_ptr<float_type> input_vector1(numElements);
    cuda_managed_ptr<float_type> input_vector2(numElements);

    // Allocate the managed output vector
    cuda_managed_ptr<float_type> output_vector(numElements);

    // a log-uniform in [1e-2, 1e3], and x / a log-uniform in [10^-1.5, 10^1.5]: this reaches every method,
    // including the tails where the non-normalised function under- or overflows.
    std::mt19937_64 gen(42);
    std::uniform_real_distribution<double> log_a(-2, 3);
    std::uniform_real_distribution<double> log_ratio(-1.5, 1.5);
    for (int i = 0; i < numElements; ++i)
    {
        double a = std::pow(10.0, log_a(gen));
        input_vector1[i] = static_cast<float_type>(a);
        input_vector2[i] = static_cast<float_type>(a * std::pow(10.0, log_ratio(gen)));
    }

    // Launch the CUDA Kernel
    int threadsPerBlock = 1024;
    int blocksPerGrid =(numElements + threadsPerBlock - 1) / threadsPerBlock;
    std::cout << "CUDA kernel launch with " << blocksPerGrid << " blocks of " << threadsPerBlock << " threads" << std::endl;

    watch w;

    cuda_test<<<blocksPerGrid, threadsPerBlock>>>(input_vector1.get(), input_vector2.get(), output_vector.get(), numElements);
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
       results.push_back(boost::math::ligamma(input_vector1[i], input_vector2[i]));
    double t = w.elapsed();
    // check the results
    for(int i = 0; i < numElements; ++i)
    {
        // The logarithm passes through zero, where relative error is meaningless, so compare absolutely there:
        float_type error = std::fabs(results[i]) < 1
            ? std::fabs(output_vector[i] - results[i]) / std::numeric_limits<float_type>::epsilon()
            : boost::math::epsilon_difference(output_vector[i], results[i]);
        if (error > 300)
        {
            std::cerr << "Result verification failed at element " << i << "!" << std::endl;
            std::cerr << "Error rate was: " << error << "eps" << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cout << "Test PASSED, normal calculation time: " << t << "s" << std::endl;
    std::cout << "Done\n";

    return 0;
}
