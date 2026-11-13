#!/usr/bin/env python3
"""Generate the pieces of a Turbulence regression case that depend on the MAPL2 capture:
the imports and tracers as data components, the root connections, and the exports to
activate. The output files are fragments, to paste into the case's YAML files
(README.md, "Rebuilding a case").

Usage: gen_case_yaml.py SPEC_RC CAPTURE_DIR OUT_DIR
  SPEC_RC      ../../Turbulence_StateSpecs.rc
  CAPTURE_DIR  the converted capture (convert_capture.sh)

Writes into OUT_DIR:
  data-for-turbulence.yaml   the unconditional imports (44), as internals
  data-tracers.yaml          the 47 TR fields of the capture, subscribed to TRDIFF
  connections_imports.yaml   root.yaml connection entries DataForTurbulence -> TURBULENCE
  history_fields.txt         the exports to activate (captured, written, not ADD2EXPORT)
  tracer_names.txt           capture name -> MAPL3 name
"""
import sys, re
import netCDF4, numpy as np

spec_path, cap_dir, out = sys.argv[1:4]

def tables(path):
    t = {}
    cat = hdr = None
    for line in open(path):
        s = line.strip()
        if s.startswith('category:'):
            cat = s.split(':')[1].strip(); hdr = None; t[cat] = []
            continue
        if s.startswith('NAME'):
            hdr = [c.strip() for c in line.split('|')]; continue
        if hdr and '|' in s and not s.startswith('#'):
            t[cat].append(dict(zip(hdr, [c.strip() for c in line.split('|')])))
    return t

T = tables(spec_path)
VDIM = {'C': 'center', 'E': 'edge', 'N': 'none'}

def field_yaml(name, std, units, dims, vloc, indent='      '):
    y = [f'{indent}{name}:',
         f'{indent}  standard_name: "{std}"',
         f'{indent}  units: "{units}"']
    if dims == 'z':
        y.append(f'{indent}  dims: "z"')
    y.append(f'{indent}  vertical_dim_spec: "{VDIM[vloc]}"')
    return '\n'.join(y)

cap = netCDF4.Dataset(f'{cap_dir}/TURBULENCE_import_before_runPhase1.nc')
capvars = set(cap.variables)

# imports: unconditional ones; all must be in the capture
imports = [r for r in T['IMPORT'] if not r['COND']]
missing = [r['NAME'] for r in imports if r['NAME'] not in capvars]
assert not missing, missing
with open(f'{out}/data-for-turbulence.yaml', 'w') as f:
    f.write('\n'.join(field_yaml(r['NAME'], r['LONG NAME'], r['UNITS'], r['DIMS'], r['VLOC'])
                      for r in imports) + '\n')
with open(f'{out}/connections_imports.yaml', 'w') as f:
    for r in imports:
        f.write(f"    - {{name: {r['NAME']}, src_comp: DataForTurbulence, dst_comp: TURBULENCE}}\n")

# tracers: TR_* of the capture, in capture order
tr = [v for v in cap.variables if v.startswith('TR_')]
names = []
with open(f'{out}/data-tracers.yaml', 'w') as f, open(f'{out}/tracer_names.txt', 'w') as g:
    for v in tr:
        name = v[3:].replace('.', '_')
        names.append(name)
        g.write(f'{v} {name}\n')
        var = cap[v]
        f.write(field_yaml(name, var.long_name, var.units, 'xyz', 'C') + '\n')
    f.write('items: [' + ', '.join(names) + ']\n')

# exports to activate: in the capture export file, not TRI_*, not internals
exp = netCDF4.Dataset(f'{cap_dir}/TURBULENCE_export_after_runPhase1.nc')
internal_names = {r['NAME'] for r in T['INTERNAL']}
export_names = {r['NAME'] for r in T['EXPORT']}
skipdims = {'Xdim', 'Ydim', 'nf', 'lev', 'edge', 'time', 'lons', 'lats', 'corner_lons',
            'corner_lats', 'contacts', 'orientation', 'anchor', 'cubed_sphere', 'ncontact'}
act, other = [], []
for v in exp.variables:
    if v in skipdims or v in exp.dimensions:
        continue
    if v.startswith('TRI_') or v in internal_names:
        other.append(v); continue
    assert v in export_names, v
    act.append(v)
with open(f'{out}/history_fields.txt', 'w') as f:
    f.write('\n'.join(act) + '\n')
print('imports', len(imports), 'tracers', len(names), 'exports to activate', len(act),
      'skipped (TRI_/internal)', len(other))
