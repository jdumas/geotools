#ifndef __COMMON_H
#define __COMMON_H

#ifdef USE_DOUBLE

// -----------------------------------------------------------------------------
#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
#error "Double precision floating point not supported by OpenCL implementation."
#endif
// -----------------------------------------------------------------------------

typedef double  real;
typedef double2 real2;
typedef double3 real3;
#define to_real3 convert_double3

#else

typedef float  real;
typedef float2 real2;
typedef float3 real3;
#define to_real3 convert_float3

#endif

// -----------------------------------------------------------------------------

typedef uint index_t;

#endif // __COMMON_H
