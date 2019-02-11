# Barcode
Bayesian Reconstruction of COsmic DEnsity fields

This repository contains both Barcode and a set of supplementary analysis tools.

[![Build Status](https://travis-ci.org/egpbos/barcode.svg?branch=master)](https://travis-ci.org/egpbos/barcode)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/db7aa18754fa4720a2d80ab47ed85e3b)](https://www.codacy.com/app/egpbos/barcode?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=egpbos/barcode&amp;utm_campaign=Badge_Grade)
[![Coverage Status](https://coveralls.io/repos/github/egpbos/barcode/badge.svg?branch=master)](https://coveralls.io/github/egpbos/barcode?branch=master)

## Citing

If you use this software, please cite it as:

Bos E. G. P., Kitaura F.-S., van de Weygaert R., 2018, preprint ([arXiv:1810.05189](https://arxiv.org/abs/1810.05189))

In bibtex format:

```bibtex
@article{Bos+2018,
      author         = "Bos, E. G. Patrick and Kitaura, Francisco-Shu and van de
                        Weygaert, Rien",
      title          = "{Bayesian cosmic density field inference from redshift
                        space dark matter maps}",
      year           = "2018",
      eprint         = "1810.05189",
      archivePrefix  = "arXiv",
      primaryClass   = "astro-ph.CO",
      journal        = {ArXiv e-prints}
}
```

Unique identifiers for citing the software itself (preferably in addition to citing the paper above) are provided through Zenodo (a unique DOI for each Barcode release) and the Astrophysics Source Code Library (all Barcode versions).

[![DOI](https://zenodo.org/badge/152633059.svg)](https://zenodo.org/badge/latestdoi/152633059) [![Join the chat at https://gitter.im/barcode_cosmo/community](https://badges.gitter.im/barcode_cosmo/community.svg)](https://gitter.im/barcode_cosmo/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
<a href="http://ascl.net/1810.002"><img src="https://img.shields.io/badge/ascl-1810.002-blue.svg?colorB=262255" alt="ascl:1810.002" /></a>

This code was previously described and applied in [conference proceedings](https://arxiv.org/abs/1611.01220) and in [Patrick's PhD thesis](https://www.rug.nl/research/portal/en/publications/clusters-voids-and-reconstructions-of-the-cosmic-web(0f7c3d17-9661-4b9f-a27c-dfac2990b844).html).

## Build

Clone the repository, `cd` into the cloned directory (`$BARCODE` below) and run the following to configure and build:

```sh
mkdir cmake-build
cd cmake-build
cmake ..
make
```

This will create binaries for barcode and the supplementary tools in the `cmake-build` directory.


## Run

Barcode must be run in the same directory as the `input.par` file.
Edit this file to change input parameters.
Then simply run with:

```
cmake-build/barcode [restart_iteration]
```

Optionally add the `restart_iteration` number when doing a restart run from existing output files.


## Development and contributing
This is an early release. Unit tests and other test codes will be added later (as mentioned in some of the code comments). Documentation as well.

Contributions are very welcome! Please don't hesitate to propose ideas for improvements in the GitHub issues or in a PR.


## License
The original contributions made as part of this code are distributed under the MIT license (see `LICENSE` file).

When compiled, this code must link to FFTW 3 and the GNU Scientific Library (GSL).
FFTW is distributed under the GNU Public License v2 or a later version, GSL under GPL v3.
This means that any redistribution of Barcode in binary form is subject to GPL v3 terms.

<!-- Barcode also depends on ncurses, which is distributed under the X11 license, a permissive license similar to the MIT license. -->
