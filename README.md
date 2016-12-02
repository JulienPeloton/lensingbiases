Lensingbiases
==

Computation of N0 and N1 biases in the context of CMB weak lensing

License: GNU License (see the LICENSE file for details) covers all files
in the lensingbiases repository unless stated otherwise.

### Before starting
The Python code has the following dependencies:
     * numpy, pylab, argparse
     * f2py
In most machines, those packages are automatically installed.
Make sure you update your PYTHONPATH to use the code.
Just add in your bashrc:
```bash
BIASPATH=/path/to/the/package
export PYTHONPATH=$PYTHONPATH:$BIASPATH
```

### Compilation of Fortran
The package is written in python, but internal operations are
done in Fortran. Therefore you need to compile it before using
the routines. We provide a setup file. Just run:
```bash
python setup.py build_ext --inplace --fcompiler=gfortran
```
Depending on system, you may need to specify a different fortran compiler.
For the puriste, we also provide a Makefile for a direct compilation.
For extended multipole range, computation can be long.
So make sure that you are using several processors by adding in your bashrc:
```bash
export OMP_NUM_THREADS=n
```

### Usage
User has to provide range of multipole, beam, noise level, etc.
In order to ease the start, we provide both a test script and a launcher:
    - LensingBiases.py
    - launch.sh
For usage on clusters, we also provide a batch (SLURM).