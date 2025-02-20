import cffi  # type: ignore

TMPFILEBASE = "pyMLINC_interface_py"

ffi = cffi.FFI()

source = """
from {} import ffi
from pyMLINC.core import pyMLINC_init, pyMLINC_run # <-- User code starts here
import traceback

@ffi.def_extern()
def pyMLINC_interface_init_py() -> int:
    try:
        # Calling out off the bridge into the python
        pyMLINC_init()
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0

@ffi.def_extern()
def pyMLINC_interface_run_py(options, in_buffer, out_buffer) -> int:
    try:
        # Calling out off the bridge into the python
        pyMLINC_run(options, in_buffer, out_buffer)
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0
""".format(TMPFILEBASE)

with open("interface.h") as f:
    data = "".join([line for line in f if not line.startswith("#")])
    data = data.replace("CFFI_DLLEXPORT", "")
    ffi.embedding_api(data)

ffi.set_source(TMPFILEBASE, '#include "interface.h"')

ffi.embedding_init_code(source)
ffi.compile(target="lib" + TMPFILEBASE + ".so", verbose=True)
