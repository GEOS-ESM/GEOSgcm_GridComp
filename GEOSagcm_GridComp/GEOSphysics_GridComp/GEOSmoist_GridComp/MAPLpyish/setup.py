#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""MAPLpyish is experimental, you shouldn't be here."""

from setuptools import find_namespace_packages, setup
from warnings import warn


warn(
    "You are installing MAPLpyish. "
    "It sounds like a terrible idea, please re-think your life choices and "
    "contact florian.g.deconinck@nasa.gov."
)

with open("README.md", encoding="utf-8") as readme_file:
    readme = readme_file.read()

setup(
    author="NASA",
    author_email="florian.g.deconinck@nasa.gov",
    python_requires=">=3.11",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Apache 2 License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.11",
    ],
    description=("MAPLpyish is experimental, you shouldn't be here."),
    install_requires=["cffi"],
    extras_require={},
    license="BSD license",
    long_description=readme,
    include_package_data=True,
    name="MAPLpyish",
    packages=find_namespace_packages(include=["MAPLpyish", "MAPLpyish.*"]),
    setup_requires=[],
    test_suite="",
    url="https://github.com/GEOS-ESM/GEOSgcm_GridComp",
    version="0.0.0",
    zip_safe=False,
)
