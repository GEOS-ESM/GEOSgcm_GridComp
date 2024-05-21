#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""pyMoist is the NDSL version of GEOS Moist physics."""

from setuptools import find_namespace_packages, setup


with open("README.md", encoding="utf-8") as readme_file:
    readme = readme_file.read()

requirements = []

test_requirements = ["pytest", "pytest-subtests", "serialbox", "coverage"]
ndsl_requirements = ["ndsl @ git+https://github.com/NOAA-GFDL/NDSL.git@2024.04.00"]
develop_requirements = test_requirements + ndsl_requirements + ["pre-commit"]

extras_requires = {
    "test": test_requirements,
    "ndsl": ndsl_requirements,
    "develop": develop_requirements,
}

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
    description="pyMoist is the NDSL version of NASA GMAO's GEOS Moist physics.",
    install_requires=requirements,
    extras_require=extras_requires,
    license="BSD license",
    long_description=readme,
    include_package_data=True,
    name="pyMoist",
    packages=find_namespace_packages(include=["pyMoist", "pyMoist.*"]),
    setup_requires=[],
    test_suite="tests",
    url="https://github.com/NOAA-GFDL/pyFV3",
    version="0.0.0",
    zip_safe=False,
)
