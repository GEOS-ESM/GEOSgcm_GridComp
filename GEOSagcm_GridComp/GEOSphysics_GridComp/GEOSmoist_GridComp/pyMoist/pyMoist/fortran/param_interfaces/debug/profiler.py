import os

import numpy as np
import yaml
from MAPL_PythonBridge import UserCode, get_MAPLPy
from MAPL_PythonBridge.types import CVoidPointer

from pyMoist.fortran.cuda_profiler import CUDAProfiler


try:
    import pyitt
except ModuleNotFoundError:
    pyitt = None


class PyProfiler(UserCode):
    def __init__(self, name: str) -> None:
        config_filename = os.getenv("GEOS_PYPROFILER_CONFIG", "")
        if os.path.exists(config_filename):
            with open(config_filename) as file:
                self._config = yaml.load(file.read(), Loader=yaml.SafeLoader)
        else:
            self._config = None

        if pyitt:
            self._pyitt_domain = pyitt.domain("GEOS Numerics")

    def _get_task_name(self, mapl_state) -> str | None:
        if not self._config:
            return None

        maplpy = get_MAPLPy()
        taskid = maplpy.get_resource("PYPROFILER_TASKID", mapl_state, default=np.int32(-1))

        if taskid not in self._config:
            print(f"pyProfiler TaskID {taskid} not found in config: {self._config}")
            return None

        return self._config[taskid]

    def init(self, mapl_state, import_state, export_state) -> None:
        taskname = self._get_task_name(mapl_state)
        if taskname is None:
            return

        if pyitt:
            pyitt.task(taskname, self._pyitt_domain).begin()

        CUDAProfiler.range_push(taskname)

    def run(self, mapl_state, import_state, export_state) -> None:
        raise RuntimeError("ITT Profiler needs to be called on .init/.finalize")

    def run_with_internal(
        self,
        mapl_state: CVoidPointer,
        import_state: CVoidPointer,
        export_state: CVoidPointer,
        internal_state: CVoidPointer,
    ) -> None:
        raise RuntimeError("ITT Profiler needs to be called on ..init/.finalize")

    def finalize(self, mapl_state, import_state, export_state) -> None:
        taskname = self._get_task_name(mapl_state)
        if taskname is None:
            return

        if pyitt:
            pyitt.task(taskname, self._pyitt_domain).end()

        CUDAProfiler.range_pop()


CODE = PyProfiler("pyProfiler")
