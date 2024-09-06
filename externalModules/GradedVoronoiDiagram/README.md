# Adaptive Graph Generation

Easily create Voronoi and Delaunay graphs from density-adaptive samples in 2D and 3D.

## Dependencies

Everything necessary can be installed via `pip install -r requirements.txt`.

**Actually required** are 
* [NumPy](https://github.com/numpy/numpy), [SciPy](https://github.com/scipy/scipy) and [drbutil](https://github.com/dbukenberger/drbutil).

**Optionally**, 
* [Mayavi](https://github.com/enthought/mayavi) for result visualization, [libigl](https://github.com/libigl/libigl-python-bindings) and [embreeX](https://github.com/trimesh/embreex) for speedup.

## Usage
Clone the repo and run `python AdaptiveGraphGenerator.py path/to/inputFile.<obj,msh>`.
Inputs should be `.obj` (2D tri) or `.msh` (3D tet) files, respectively.
Scalar density values for the vertices are expected in a separate `inputFile.scl` file.
See `/data` for examples.

Additional command line parameters are:
```
--s	Scale factor for the Poisson disk radii (multiplied with the bounding box diagonal length, default is 0.05).
--r	Ratio of smallest to largest Poisson radius (default is 0.25).
--P	Show Poisson samples (requires mayavi).
--G	Show graphs (requires mayavi).
```

Generated outputs are named like the input files.
Filenames further give the number of samples `n`, as well as used parameters `s` and `r`.