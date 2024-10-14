#!/usr/bin/env python

import os
import subprocess
import pathlib

try:
    subprocess.run(["ml", "SMTStack/2024.04.00"])
except FileNotFoundError:
    print("Could not load the module, trusting you know what you are doing")

this_dir_path = os.path.dirname(os.path.realpath(__file__))
moist_comp = os.path.abspath(os.path.join(this_dir_path, "../../../"))
print(f"Running all .SER in {moist_comp}")

ppser_args = ["--verbose", "--ignore-identical", "-m", "utils_ppser_kbuff"]
ppser_script = os.path.abspath(f"{os.getenv('SERIALBOX_ROOT')}/python/pp_ser/pp_ser.py")
print(f"Executing ppser from {ppser_script}")


for path in pathlib.Path(moist_comp).rglob("*.SER"):
    cmd = ["python", ppser_script]
    cmd.extend(ppser_args)
    cmd.extend(["-o", str(path.absolute())[:-4], str(path.absolute())])
    subprocess.run(cmd)
