#!/usr/bin/env python3
"""Copy a converted (new cubed-sphere) capture file, keeping the grid variables
and the listed data variables, optionally renamed.

Usage: subset_nc.py IN OUT name[=newname] ...   (or @FILE with one 'old new' per line)
       subset_nc.py IN OUT --all-but-prefix TR_
"""
import sys
import netCDF4

src, dst, *args = sys.argv[1:]
inp = netCDF4.Dataset(src)
grid = {'Xdim', 'Ydim', 'nf', 'lev', 'edge', 'time', 'lons', 'lats', 'corner_lons',
        'corner_lats', 'contacts', 'orientation', 'anchor', 'cubed_sphere'}
data = [v for v in inp.variables if v not in grid and v not in inp.dimensions]

pairs = []
if args and args[0] == '--all-but-prefix':
    pairs = [(v, v) for v in data if not v.startswith(args[1])]
else:
    for a in args:
        if a.startswith('@'):
            pairs += [tuple(l.split()) for l in open(a[1:]) if l.strip()]
        else:
            o, _, n = a.partition('=')
            pairs.append((o, n or o))
missing = [o for o, _ in pairs if o not in inp.variables]
assert not missing, missing

out = netCDF4.Dataset(dst, 'w', format='NETCDF4')
out.setncatts({k: inp.getncattr(k) for k in inp.ncattrs()})
for name, dim in inp.dimensions.items():
    out.createDimension(name, None if dim.isunlimited() else len(dim))

def copy(old, new):
    v = inp.variables[old]
    fill = v.getncattr('_FillValue') if '_FillValue' in v.ncattrs() else None
    w = out.createVariable(new, v.datatype, v.dimensions, fill_value=fill)
    w.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != '_FillValue'})
    w[...] = v[...]

for g in inp.variables:
    if g in grid:
        copy(g, g)
for o, n in pairs:
    copy(o, n)
out.close()
print(f'{dst}: {len(pairs)} variables')
