-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
Crystalline Ice Finder v0.5
-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
USAGE: find_xtal_ice.py --i <path to directory containing micrographs/particle stacks or image mrc file> --apix <angstroms per pixel>
If the input images are stacks add '--stack' to the input line

If entering a directory make sure to include the trailing / on the directory name; the program adds *.mrc* to the end of it to find both .mrc or .mrcs files

The script operates on .mrc single images (micrographs) or .mrcs stacks (relion particle stacks) but not both at the same time

Required Modules:
glob
sys
os
matplotlib.pyplot
scipy fftpack
scipy misc
numpy
struct

The program may throw a warning:

/fbs/emsoftware2/env-modules/install/anaconda/2018.12/lib/python2.7/site-packages/scipy/fftpack/basic.py:160: FutureWarning: Using a non-tuple sequence for multidimensional indexing is deprecated; use `arr[tuple(seq)]` instead of `arr[seq]`. In the future this will be interpreted as an array index, `arr[np.array(seq)]`, which will result either in an error or a different result.
  z[index] = x

This is due to using an older version scipy - update you install of scipy
