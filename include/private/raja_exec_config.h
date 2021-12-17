#pragma once

/**
 * @brief Execution policy logic for raja-based code
 * @file include/exec_policy_config.hpp
 * @author Asher Mancinelli <asher.mancinelli@pnnl.gov>
 *
 * @todo Our CUDA and HIP parallel execution policies currently use 128 threads
 * no matter the situation. In the furutre, this should be parameterized via
 * CMake.
 */

#include <exago_config.h>

#if !defined(EXAGO_ENABLE_RAJA)
#error "`raja_exec_config.hpp` was included, but RAJA is not enabled!"
#endif

#include <RAJA/RAJA.hpp>

#if defined(EXAGO_ENABLE_GPU)

/* Annote all lambdas with `__device__` qualifier when running on device */
#define RAJA_LAMBDA [=] __device__

/* Determine whether CUDA, HIP, or OpenMP policies should be used */
#if defined(EXAGO_ENABLE_CUDA)

using exago_raja_exec = RAJA::cuda_exec<128>;
using exago_raja_reduce = RAJA::cuda_reduce;
using exago_raja_atomic = RAJA::cuda_atomic;

#elif defined(EXAGO_ENABLE_HIP)

using exago_raja_exec = RAJA::hip_exec<128>;
using exago_raja_reduce = RAJA::hip_reduce;
using exago_raja_atomic = RAJA::hip_atomic;

#endif // defined(EXAGO_ENABLE_CUDA)

#else // EXAGO_ENABLE_GPU

/* If not running on the device, just use by-value capture cluase by default */
#define RAJA_LAMBDA [=]

using exago_raja_exec = RAJA::omp_parallel_for_exec;
using exago_raja_reduce = RAJA::omp_reduce;
using exago_raja_atomic = RAJA::omp_atomic;

#endif // EXAGO_ENABLE_GPU
