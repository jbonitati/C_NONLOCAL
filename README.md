# Multichannel R-Matrix Solver
[git-repo-url]
> Authors: Joey Bonitati, Weichuan Li, Gregory Potel, Filomena Nunes

Known Issues:
- The normalization of the channels in the coupled channel case mysteriously depends on the number of basis functions being used

## Synopsis

This program uses the calculable R-matrix method described in [Descouvement-2010] to compute wave functions for multichannel collisions involving local and non-local optical potentials.

## Installation

This program requires the installation of the [Armadillo], [Boost], and [GSL] libraries. If these libraries are not installed to automatically linked location, you can add -I flags to the LDFLAGS variable in the makefile.

To compile use the command 
```
$ make
```
And to execute the program use the command
```
$ ./run
```

## Tests

Example input files can be found in the "examples" directory, and outputs are in the "output" directory.

The program "examples/makeini.py" can be used to generate a sample input file for an arbitrary number of channels. For example:
```
$ python examples/makeini.py config.ini 5
```
will overwrite the "config.ini" file with input parameters for a 5 channel problem, which can then be modified in a text editor.

[//]: #
   [git-repo-url]: <https://github.com/jbonitati/C_NONLOCAL.git>
   [Descouvement-2010]: <https://arxiv.org/abs/1001.0678>
   [Armadillo]: <http://arma.sourceforge.net/>
   [Boost]: <http://www.boost.org/>
   [GSL]: <https://www.gnu.org/software/gsl/>
/README.md>
