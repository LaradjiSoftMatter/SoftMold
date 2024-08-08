#include <cuda.h>
#include <curand_kernel.h>

#ifndef MPD_ERRORS
#define MPD_ERRORS

#define CUDA_Kernel_Errors() {\
		cudaError_t CUDA_RT_ERROR = cudaGetLastError();\
		if(CUDA_RT_ERROR!=cudaSuccess) {\
                std::string error=cudaGetErrorName(CUDA_RT_ERROR);\
                error+=" on line ";\
                error+=std::to_string(__LINE__);\
                error+=" in file ";\
                error+=__FILE__;\
                error+=" in function ";\
                error+=__func__;\
                throw std::runtime_error(error);}}

#define CUDA_API_Errors(CUDA_FUNCTION) {\
		cudaError_t CUDA_RT_ERROR=(CUDA_FUNCTION);\
		if(CUDA_RT_ERROR!=cudaSuccess) {\
                std::string error=cudaGetErrorName(CUDA_RT_ERROR);\
                error+=" on line ";\
                error+=std::to_string(__LINE__);\
                error+=" in file ";\
                error+=__FILE__;\
                error+=" in function ";\
                error+=__func__;\
		error+=" from CUDA API ";\
                throw std::runtime_error(error);}}

#define CUDA_API_Warnings(CUDA_FUNCTION) {\
		cudaError_t CUDA_RT_ERROR=(CUDA_FUNCTION);\
		if(CUDA_RT_ERROR!=cudaSuccess) {\
                std::string error=cudaGetErrorName(CUDA_RT_ERROR);\
                error+=" on line ";\
                error+=std::to_string(__LINE__);\
                error+=" in file ";\
                error+=__FILE__;\
                error+=" in function ";\
                error+=__func__;\
		error+=" from CUDA API ";\
                std::cerr << error << std::endl;}}

static const char *curandErrorString(curandStatus_t error)
{
    switch (error)
    {
        case CURAND_STATUS_SUCCESS:
            return "CURAND_STATUS_SUCCESS";

        case CURAND_STATUS_VERSION_MISMATCH:
            return "CURAND_STATUS_VERSION_MISMATCH";

        case CURAND_STATUS_NOT_INITIALIZED:
            return "CURAND_STATUS_NOT_INITIALIZED";

        case CURAND_STATUS_ALLOCATION_FAILED:
            return "CURAND_STATUS_ALLOCATION_FAILED";

        case CURAND_STATUS_TYPE_ERROR:
            return "CURAND_STATUS_TYPE_ERROR";

        case CURAND_STATUS_OUT_OF_RANGE:
            return "CURAND_STATUS_OUT_OF_RANGE";

        case CURAND_STATUS_LENGTH_NOT_MULTIPLE:
            return "CURAND_STATUS_LENGTH_NOT_MULTIPLE";

        case CURAND_STATUS_DOUBLE_PRECISION_REQUIRED:
            return "CURAND_STATUS_DOUBLE_PRECISION_REQUIRED";

        case CURAND_STATUS_LAUNCH_FAILURE:
            return "CURAND_STATUS_LAUNCH_FAILURE";

        case CURAND_STATUS_PREEXISTING_FAILURE:
            return "CURAND_STATUS_PREEXISTING_FAILURE";

        case CURAND_STATUS_INITIALIZATION_FAILED:
            return "CURAND_STATUS_INITIALIZATION_FAILED";

        case CURAND_STATUS_ARCH_MISMATCH:
            return "CURAND_STATUS_ARCH_MISMATCH";

        case CURAND_STATUS_INTERNAL_ERROR:
            return "CURAND_STATUS_INTERNAL_ERROR";
    }

    return "<unknown>";
}


#define CURAND_Errors(CURAND_FUNCTION) do { \
		curandStatus_t CURAND_ST_ERROR=(CURAND_FUNCTION); \
		if(CURAND_ST_ERROR != CURAND_STATUS_SUCCESS) { \
        	std::string error=curandErrorString(CURAND_ST_ERROR);\
                error+=" on line ";\
                error+=std::to_string(__LINE__);\
                error+=" in file ";\
                error+=__FILE__;\
                error+=" in function ";\
                error+=__func__;\
		error+=" from CURAND API ";\
		throw std::runtime_error(error);}} while(0)

#endif
