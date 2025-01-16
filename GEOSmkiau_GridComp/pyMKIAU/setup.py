#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""pyMKIAU - python sub-component of GEOS MKIAU."""

from setuptools import find_namespace_packages, setup


with open("README.md", encoding="utf-8") as readme_file:
    readme = readme_file.read()

setup(
    author="NASA",
    python_requires=">=3.11",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Apache 2 License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.11",
    ],
    description=("pyMKIAU - python sub-component of GEOS MKIAU."),
    install_requires=[],
    extras_require={},
    long_description=readme,
    include_package_data=True,
    name="pyMKIAU",
    packages=find_namespace_packages(include=["pyMKIAU", "pyMKIAU.*"]),
    setup_requires=[],
    url="https://github.com/GEOS-ESM/GEOSgcm_GridComp",
    version="0.0.0",
    zip_safe=False,
)
