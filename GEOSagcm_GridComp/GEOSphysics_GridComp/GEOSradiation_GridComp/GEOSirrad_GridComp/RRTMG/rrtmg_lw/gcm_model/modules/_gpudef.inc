#ifdef _CUDA
#define _gpudev  ,device
#define _gpudeva ,device,allocatable
#define _gpudevanp ,device,allocatable
#define _gpucon  ,constant
#define _gpuker  attributes(global)
#define _gpuked  attributes(device)
#define _gpuchv  <<<dimGrid,dimBlock>>>
#define _cpus    
#define _cpusnp
#else
#define _gpudev
#define _gpudeva ,pointer
#define _gpudevanp ,allocatable
#define _gpucon
#define _gpuker
#define _gpuked
#define _gpuchv
#define _cpus    ,target
#define _cpusnp  
#endif