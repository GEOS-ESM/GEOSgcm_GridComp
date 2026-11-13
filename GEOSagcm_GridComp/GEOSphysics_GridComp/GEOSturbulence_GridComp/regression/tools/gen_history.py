#!/usr/bin/env python3
"""Writes history.yaml of a case's run1 or run2 from the export list that
gen_case_yaml.py wrote (history_fields.txt).

Usage: gen_history.py CASE_RUN_DIR RUN FIELDS_FILE [EXCLUDE ...]
  CASE_RUN_DIR   e.g. .../regression/turbulence-default-sa/run1 (history.yaml is written there)
  RUN            run1 or run2 (run1 also drops the 9 exports only Run2 writes)
  EXCLUDE        exports this configuration never writes, or that must not be requested
                 (README.md, "Which exports to request")
"""
import sys

rundir, run, ff, *excl = sys.argv[1:]
RUN2_ONLY = ['DPDTTRB', 'INTDIS', 'QT', 'SAFUPDATE', 'SL', 'SRFDIS', 'TOPDIS', 'WQT', 'WSL']
names = [l.strip() for l in open(ff) if l.strip()]
names = [n for n in names if n not in excl and not (run == 'run1' and n in RUN2_ONLY)]

hdr = f"""# Activates (requests) the TURBULENCE exports the MAPL2 run had allocated
# (requested by History or by other components, or ALLOC'd): those in the MAPL2
# export capture, less the TRI bundle members and the fields that are internals
# in MAPL3 (always active), and less the ones that are not requested in this
# configuration: {', '.join(sorted(excl)) or 'none'} (MAPL3 would leave a
# requested export that is never written uninitialized, MAPL2 zeroed it).
"""
if run == 'run1':
    hdr += """# Run1 only: also less the 9 that only Run2 fetches and writes (DPDTTRB, INTDIS,
# QT, SAFUPDATE, SL, SRFDIS, TOPDIS, WQT, WSL; MAPL2's capture after Run1 holds
# the previous step's values).
"""
hdr += """# Requesting a field changes what MAPL2's ALLOC= and associated() tests do, so
# this list reproduces the MAPL2 run. A collection that never writes (start
# after the cap stop) still connects its fields. No geom/vertical_grid: native
# grid.
mapl:
  setServices:
    sharedObj: libMAPL.history

active_collections:
  - turbulence

time_specs:
  never: &never
    frequency: PT20M
    start: "3000-01-01T00:00:00"

collections:
  turbulence:
    template: "%c.%y4%m2%d2_%h2%n2z.nc4"
    time_spec: *never
    var_list:
"""
w = max(len(n) for n in names) + 2
body = ''.join(f"      {n + ':':<{w}} {{source: TURBULENCE.{n}}}\n" for n in sorted(names))
with open(f'{rundir}/history.yaml', 'w') as f:
    f.write(hdr + body)
print(f'{rundir}/history.yaml: {len(names)} exports')
