# Fortran - Python bridge prototype

Nomenclatura: we call the brige "fpy" and "c", "f" and "py" denotes functions in their respective language.

Building: you have to pass `-DBUILD_PYMKIAU_INTERFACE=ON` to your `cmake` command to turn on the interface build and execution.

## Pipeline

Here's a quick rundown of how a buffer travels through the interface and back.
  
- From Fortran in `GEOS_MKIAUGridComp:488` we call `pyMKIAU_interface_f_run` with the buffer passed as argument
- This pings the interface, located at `pyMKIAU/interface/interface.f90`. This interface uses the `iso_c_binding` to marshall the parameters downward (careful about the user type, look at the code)
- Fortran then call into C at `pyMKIAU/interface/interface.c`. Those functions now expect that a few `extern` hooks have been made available on the python side, they are define in `pyMKIAU/interface/interface.h`
- At runtime, the hooks are found and code carries to the python thanks to cffi. The .so that exposes the hooks is in `pyMKIAU/interface/interface.py`. Within this code, we: expose extern functions via `ffi.extern`, build a shared library to link for runtime and pass the code down to the `pyMKIAU` python package which lives at `pyMKIAU/pyMKIAU`
- In the package, the `serservices` or `run` function is called.

## Fortran <--> C: iso_c_binding

We leverage Fortan `iso_c_binding` extension to do conform Fortran and C calling structure. Which comes with a bunch of easy type casting and some pretty steep potholes.
The two big ones are:

- strings need to be send/received as a buffer plus a length,
- pointers/buffers are _not_ able to be pushed into a user type.

## C <->Python: CFFI based glue

The interface is based on CFFI which is reponsible for the heavy lifting of

- spinning a python interpreter
- passing memory between C and Python without a copy

## Running python

The last trick is to make sure your package is callable by the `interface.py`. Basically your code has to be accessible by the interpreter, be via virtual env, conda env or PYTHONPATH. The easy way to know is that you need to be able to get into your environment and run in a python terminal

```python
from pyMKIAU.core import pyMKIAU_init
pyMKIAU_init()
```
