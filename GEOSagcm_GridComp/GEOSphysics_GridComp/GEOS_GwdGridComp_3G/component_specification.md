Specification of a _leaf_ gridded component involves the following actions:

- Extracting parameter values from a configuration
  - YAML (preferred)
  - ESMF (deprecated)
  - Namelist (strongly discouraged)
  
- Specification of field attributes for import, export, and internal state
  - `#include` from ACG (preferred)
  - Explicit `add*Spec()` (deprecated)
  
- Specification of entry points:  init/run/finalize/... with optional phases

- Specification of grid
  - inherit from parent (usual case)
  - define local grid

Parent (non-leaf) components have the actions above as well as the
following additional configuration actions:

- Specification of child components
  - YAML
  
- Specification of connections among children
  - YAML
  
To support conditional connections, the same YAML file should be used
for the entire specification.

Must still have a mechanism for user procedure in case config file
alone cannot handle logic.
