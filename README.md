# barcode
Bayesian Reconstruction of COsmic DEnsity fields

This repository contains both Barcode and a set of supplementary analysis tools.


# Build

Clone the repository, `cd` into the cloned directory (`$BARCODE` below) and run the following to configure and build:

```sh
mkdir cmake-build
cd cmake-build
cmake ..
make
```

This will create binaries for barcode and the supplementary tools in the `bin` directory.

# Run

Barcode must be run in the same directory as the `input.par` file.
Edit this file to change input parameters.
Then simply run with:

```
$BARCODE/bin/barcode [restart_iteration]
```

Optionally add the `restart_iteration` number when doing a restart run from existing output files.


# MacOS
The MacOS version has no OpenMP support. Configure the project with `cmake` options `-DMULTITHREAD=OFF -DMULTITHREAD_FFTW=OFF` to build on macOS.


# Development and contributing
This is an early release. Unit tests and other test codes will be added later (as mentioned in some of the code comments). Documentation as well.

Contributions are very welcome! Please don't hesitate to propose ideas for improvements in the GitHub issues or in a PR.


# License
The original contributions made as part of this code are distributed under the MIT license (see `LICENSE` file).

Barcode makes use of an early version of the [Planck LevelS toolkit](https://sourceforge.net/projects/planck-ls/), which is distributed under the GNU Public License v2 (see `planck/LICENSE` file).

This code depends on FFTW 3, the GNU Scientific Library, both distributed under the GPL as well.

These three dependencies mean that any redistribution of Barcode in binary form is subject to GPLv2 terms.

Barcode also depends on ncurses, which is distributed under the X11 license, a permissive license similar to the MIT license.
