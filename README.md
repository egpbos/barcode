# Barcode
Bayesian Reconstruction of COsmic DEnsity fields

This repository contains both Barcode and a set of supplementary analysis tools.

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

A Zenodo DOI is provided for each release to cite the software itself (preferably in addition to citing the paper above). Latest release:

[![DOI](https://zenodo.org/badge/152633059.svg)](https://zenodo.org/badge/latestdoi/152633059)

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

Barcode makes use of an early version of the [Planck LevelS toolkit](https://sourceforge.net/projects/planck-ls/), which is distributed under the GNU Public License v2 (see `planck/LICENSE` file).

This code depends on FFTW 3, the GNU Scientific Library, both distributed under the GPL as well.

These three dependencies mean that any redistribution of Barcode in binary form is subject to GPLv2 terms.

Barcode also depends on ncurses, which is distributed under the X11 license, a permissive license similar to the MIT license.
