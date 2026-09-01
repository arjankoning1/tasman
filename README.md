# TASMAN

TASMAN is a statistical software package for the TALYS nuclear model code. The most important features are:

- uncertainty distributions and covariance matrices for TALYS results, through Monte Carlo sampling of TALYS nuclear model parameters
- automatic optimization ('search') of TALYS nuclear model parameters to experimental nuclear reaction data and data from nuclear data libraries
- Total Monte Carlo: generation of a statistical ensemble of complete nuclear data libraries for uncertainty propagation
- parameter sensitivity profiles for TALYS output

## Documentation and reference

A description of the code and its options can be found in the [TASMAN Tutorial (pdf)](https://github.com/arjankoning1/tasman/blob/main/doc/tasman.pdf).

The reference to be used for TASMAN is:

A.J. Koning, D. Rochman, J.-Ch. Sublet, N. Dzysiuk, M. Fleming, and S. van der Marck, *TENDL: Complete Nuclear Data Library for innovative Nuclear Science and Technology*, Nuclear Data Sheets 155, 1 (2019).

## Installation

### Prerequisites

The following are the prerequisites for compiling TASMAN:

- git (only if the package is downloaded via GitHub)
- GNU make
- a recent Fortran compiler, such as GNU Fortran (gfortran)
- a successful installation of the TALYS nuclear model code
- for Total Monte Carlo: a successful installation of TEFAL for ENDF-6 formatting

### Downloads

To download TASMAN, you can use one of the following options.

#### 1. Download the entire tar file (frozen version TASMAN-2.2)

This is available at the the [TALYS page](https://nds.iaea.org/talys/), and can be retrieved by clicking on the download link or
```bash
curl -LO https://nds.iaea.org/talys/codes/tasman.tar
tar zxf tasman.tar
```

#### 2. Using git (latest beta version)

```bash
git clone https://github.com/arjankoning1/tasman.git
```

### Additional data

Full use of TASMAN, including optimization to experimental or evaluated nuclear data, requires additional databases.

For optimization to experimental data, install EXFORTABLES:

```bash
curl -LO https://nds.iaea.org/talys/exfortables.tar
tar zxf exfortables.tar
```

For optimization to evaluated nuclear data, download the required nuclear data library archives from:

```text
https://nds.iaea.org/talys/
```

TASMAN derives these external paths from the parent directory of `TASMAN_DIR`. A typical installation layout is therefore:

```text
parent_directory/
├── talys/
├── tefal/
├── tasman/
├── exfortables/
└── libraries/
```

The optional photon strength function database is expected as:

```text
parent_directory/PSF/Photo/
```

TASMAN does not use a shared parent-level `bin/` directory. External executables such as `talys`, `tefal`, `tares`, `tafis` and `tanes` are resolved through the user's `PATH`. Therefore, any required executable must be available through `PATH`, for example:

```bash
export PATH="/path/to/talys/bin:/path/to/tefal/bin:$PATH"
```

### Installation instructions

#### 1. For the tar file (frozen version TASMAN-2.2)

```bash
cd tasman
./install_tasman.bash
```

An alternative option is:

```bash
cd tasman/source
make
```

The above will invoke the default compiler `gfortran`.

#### 2. For the git version (latest beta version)

```bash
cd tasman
./install_tasman.bash
```

which automatically executes the `Makefile` in `tasman/source`. At the end, `install_tasman.bash`
will print the recommended shell configuration.

An alternative option is:

```bash
cd tasman/source
make
```

For the git version, the default compiler is `gfortran`. When `gfortran` is used and no `FFLAGS` are supplied, the Makefile uses:

```text
-w -O3 -ffp-contract=off
```

For other compilers, no default compiler flags are imposed.

The compiler and compilation options can be passed to the Makefile through `install_tasman.bash`. For example:

```bash
# GNU Fortran
./install_tasman.bash FC=gfortran FFLAGS="-O3 -ffp-contract=off"

# Intel Fortran
./install_tasman.bash FC=ifx FFLAGS="-O3"
```

The above will produce a `tasman` executable in the `tasman/bin` directory.

Set `TASMAN_DIR` to the TASMAN installation directory. This variable is required unless the fallback path in `source/machine.f90` has been set manually. For example:

```bash
export TASMAN_DIR="/Users/koning/tasman"
```

If you want to run `tasman` from anywhere, add the TASMAN `bin` directory to `PATH`:

```bash
export PATH="$TASMAN_DIR/bin:$PATH"
```

To include your name in the output files, set:

```bash
export TASMAN_USER="Your Name"
```

These lines can be added to your shell configuration file, for example `~/.zshrc` or `~/.profile`.

If setting `TASMAN_DIR` is not possible on a particular system, edit `code_dir` in `source/machine.f90` and rebuild TASMAN.

## The TASMAN package

The `tasman/` directory contains the following directories and files:

- `README.md` this README file
- `LICENSE` the License file
- `install_tasman.bash` installation script
- `source/` the Fortran source code of TASMAN and the Makefile
- `bin/` the `tasman` executable after successful installation
- `misc/` files with TASMAN input settings and the EXFOR outlier table
- `parameters/` files with probability distributions for Bayesian Monte Carlo
- `doc/` the tutorial in PDF format
- `samples/` the input and output files of the sample cases, and the `verify` script used to run the sample cases

In total, about 2.5 GB of free disk space is required to install TASMAN.

### Miscellaneous options

The standard TASMAN package is sufficient for normal use. A less frequently used feature requires the additional `PSF/` database with experimental photon strength functions for fitting TALYS nuclear model parameters.

## Sample cases

The sample cases provide examples of the use of TASMAN and can be used to verify a successful installation. The `samples/` directory contains various sample cases with a subdirectory `org/` containing the reference results and a subdirectory `new/` for results produced by the user.

The TASMAN sample cases assume that TALYS and TEFAL are installed in sibling directories:

```text
parent_directory/
├── talys/
├── tefal/
└── tasman/
```

The TALYS and TEFAL executables must also be available through `PATH` for normal TASMAN runs.

To run the sample cases:

```bash
cd samples
./verify
```

From the top-level TASMAN directory, the same test can be started with:

```bash
make -C source check
```

`make check` automatically sets `TASMAN_DIR` and the sibling `TALYS_DIR` and `TEFAL_DIR`, and adds the TALYS and TEFAL
executables to `PATH` for the test.

You may create your own input file, for example `tasman.inp`, after which TASMAN works as follows:

```bash
tasman < tasman.inp > tasman.out
```

assuming that `tasman/bin` has been added to `PATH`.

## License and Copyright

This software is distributed and copyrighted according to the [LICENSE](LICENSE) file.
