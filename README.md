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

## Inputs

The input file requires all of the following values:
- Settings
-- output_file: name of the .txt file to output to
-- num_channels: number of channels to consider (you must have a separate section for each channel)
-- entrance_channel: the number (starting from 1) of the channel to use as the entrance channel
- Numerical
-- Projectile_mass_number
-- Target_mass_number
-- Projectile_proton_number
-- Target_proton_number
-- Coulomb_radius: radius (fm) used in the coulomb potential if applicable
-- Total_energy: E (MeV) in the center of mass frame
-- Basis_size: number of Basis functions to use
-- Step_size: step (fm) at which to calculate the wave functions
-- R_max: max value to calculate the wave function
-- Channel_radius: a (fm)
- Channelx (replace x with the number >= 1)
-- Angular_momentum: l
-- Spin: m
-- Total_angular_momentum: j
-- Energy: Channel energy (MeV)
- couplingxy (replace xy with two numbers >= 1)
-- beta: coupling strength
-- V: MeV
-- r: fm
-- a: fm
- local: Volume, surface, and spin-orbit parameters for optical potential (all in MeV or fm) (make V=0 for any unwanted potentials)
-- Vv, av, rv, Wv, awv, rwv, Vd, ad, rd, Wd, awd, rwd, Vso, aso, rso, Wso, awso, Rwso
- nonlocal: same as local, but with nonlocality (make beta=0 for no nonlocality)
-- beta, Vv, av, rv, Wv, awv, rwv, Vd, ad, rd, Wd, awd, rwd, Vso, aso, rso, Wso, awso, Rwso

[//]: #
   [git-repo-url]: <https://github.com/jbonitati/C_NONLOCAL.git>
   [Descouvement-2010]: <https://arxiv.org/abs/1001.0678>
   [Armadillo]: <http://arma.sourceforge.net/>
   [Boost]: <http://www.boost.org/>
   [GSL]: <https://www.gnu.org/software/gsl/>
/README.md>
