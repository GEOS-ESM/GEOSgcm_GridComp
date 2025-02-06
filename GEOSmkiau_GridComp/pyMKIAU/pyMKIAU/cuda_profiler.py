import time
from typing import Dict, List


# Conditional cupy import for non-GPU machines
try:
    import cupy as cp
except ModuleNotFoundError:
    cp = None

# Run a deviceSynchronize() to check
# that the GPU is present and ready to run
if cp is not None:
    try:
        cp.cuda.runtime.deviceSynchronize()
        GPU_AVAILABLE = True
    except cp.cuda.runtime.CUDARuntimeError:
        GPU_AVAILABLE = False
else:
    GPU_AVAILABLE = False


class CUDAProfiler:
    """Leverages NVTX & NSYS to profile CUDA kernels."""

    def __init__(self, label: str) -> None:
        self.label = label

    def __enter__(self):
        if GPU_AVAILABLE:
            cp.cuda.runtime.deviceSynchronize()
            cp.cuda.nvtx.RangePush(self.label)

    def __exit__(self, _type, _val, _traceback):
        if GPU_AVAILABLE:
            cp.cuda.runtime.deviceSynchronize()
            cp.cuda.nvtx.RangePop()

    @classmethod
    def sync_device(cls):
        if GPU_AVAILABLE:
            cp.cuda.runtime.deviceSynchronize()

    @classmethod
    def start_cuda_profiler(cls):
        if GPU_AVAILABLE:
            cp.cuda.profiler.start()

    @classmethod
    def stop_cuda_profiler(cls):
        if GPU_AVAILABLE:
            cp.cuda.profiler.stop()

    @classmethod
    def mark_cuda_profiler(cls, message: str):
        if GPU_AVAILABLE:
            cp.cuda.nvtx.Mark(message)


class TimedCUDAProfiler(CUDAProfiler):
    def __init__(self, label: str, timings: Dict[str, List[float]]) -> None:
        super().__init__(label)
        self._start_time = 0
        self._timings = timings

    def __enter__(self):
        super().__enter__()
        self._start_time = time.perf_counter()

    def __exit__(self, _type, _val, _traceback):
        super().__exit__(_type, _val, _traceback)
        t = time.perf_counter() - self._start_time
        if self.label not in self._timings:
            self._timings[self.label] = [t]
        else:
            self._timings[self.label].append(t)
