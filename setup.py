from setuptools import setup, find_packages


setup(
    name = "vgrid",
    version = "1.0.0",
    license = "BSD 2-Clause License",
    packages = find_packages(),
    python_requires=">=3.0",
    install_requires=[
        "numpy",
        "scipy",
    ],
    description="Module for incremental gridding of XYZ data.",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: GIS"
    ],
    keywords = "statistics GIS multibeam hydrography",
    author = "Val Schmidt, Eric Younkin",
    author_email="vschmidt@ccom.unh.edu, eric.g.younkin@noaa.gov"

)
