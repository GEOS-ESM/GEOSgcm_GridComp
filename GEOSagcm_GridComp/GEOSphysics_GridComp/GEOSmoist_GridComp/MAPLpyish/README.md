# MAPLpyish: Exposing MAPL to Python

This is a prototype that exposes the MAPL Fortran to python via `cffi`. It brings the state into Python as a blind pointer and read-only memory repository (unenforced) and present an API to operate on it

## Example

A small example take `n_modes` from the `AERO` state is given. See commit `56a8470201ecac828ffaae590151a29d96bafa08`

## [Dev] TODO

- [ ] Make `scalar` and `pointer` interfaces. Scalars should probably be broken down further into types. Pointer can all go through `c_ptr`. As bridle as it is to go through blind pointer, we don't have a way ATM to enforce typing throughout
- [ ] Refactor: use `cffi` API mode rather ABI to enhance safety and portability
- [ ] Freeze MAPL state into it's own class and make it read-only

## Strategy

### Fortran <-> C

See file `c_bridge_to_MAPL.F90`.

First we have to expose the Fortran to C. This is easily done via the `iso_c_binding`, but there's some restrictions to have the python working later or due to C/Fortran differences.

#### Pointers

Fortran & C don't agree on pointer addressing. Functions to pass from one definition to the others:

- `c_f_pointer` to modify a C pointer to Fortran
- `c_loc` to pass down a Fortran pointer down to C

#### String

Passing a string requires to pass the `char*` pointer and the length of the string. It's then reconstructed Fortran side with the combination:

```fortran
type(c_ptr), intent(in), value :: name_c_ptr
integer(c_int), intent(in), value :: name_len
character(len=name_len,kind=c_char), pointer :: name
...
call c_f_pointer(name_c_ptr, name)
```

#### ESMF/MAPL state

The ESMF/MAPL state is considered an unmutable black box. _Nobody but Fortran should touch it_. To make sure it's the case, we pass it as a `void*` all the way down and up.

Pass with:

```fortran
c_loc(state)
```

Recover with:

```fortran
c_f_pointer(state_c_ptr, state)
```

### Python <-> C

See file `api.py`.

We rely on `cffi` for our work. This POC uses the ABI mode, which means we need to load the DLL/SO first. This is probably not the final way to do it, since `cffi` has also an API mode that is probably a better idea (and is used for the other interfacing we do).

#### `cffi` ABI mode

- Load the .so with `ffi.dlopen`
- Defined to `ffi` the signature of the C function with `cdef` (must match binding of Fortran)
- Every call up in the API then talks directly to the fortan_c bridge lib with marshalling data

#### String

We need to pass strings as `char*` + length.

```python
name = "a string"
name_as_pchar = self.ffi.new("char[]", name.encode())
lib.fortran_func(name_as_pchar, len(name))
```

## Some refs

How to carry struct via C ptr: <https://community.intel.com/t5/Intel-Fortran-Compiler/Fortran-C-interopability-passing-void-pointers/td-p/893892>
How to string: <https://fortran-lang.discourse.group/t/passing-strings-from-and-to-fortran-dll-using-ctypes/6482/2>
How to field: <https://www.matecdev.com/posts/fortran-in-python.html>
