-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
Crystalline Ice Finder v0.5
-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
USAGE: find_xtal_ice.py --i <path to directory containing micrographs/particle stacks or image mrc file> --apix <angstroms per pixel>
If the input images are stacks add '--stack' to the input line

If entering a directory make sure to include the trailing / on the directory name; the program adds *.mrc* to the end of it to find both .mrc or .mrcs files

The script operates on .mrc single images (micrographs) or .mrcs stacks (relion particle stacks) but not both at the same time

-=-=-=-=-=-=-=-=-=-=-=-=-=
particle filtering scripts
-=-=-=-=-=-=-=-=-=-=-=-=-=
Either select particles below a certain crystalline ice threshold of lowpass filter the paticles with too much ice to 4.0A resolution.

USAGE: particles-filter-by-ice.py --parts <particle starfile> --thresh <threshold>
If the log file from xice_finder is called somthing different than xice_find.log specify it with --ice <path to logfile>

An interesting idea but the utility of this is questionable... it doesn't seem to result in any gains in resolution generally.  If ice contamination is bad enough to cause significant drops in resolution at the ice diffraction frequencys your dataset is probably shit anyway.

-=-=-=-=-=-=-=-=-
Required Modules:
-=-=-=-=-=-=-=-=-
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
