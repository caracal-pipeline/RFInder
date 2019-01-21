# RFInder


This is a set of tools that have been developed in preparation of the [Apertif surveys](
https://www.astron.nl/astronomy-group/apertif/science-projects/apertif-science-projects).

The main function of `rfinder` is to identify the presence of RFI in an observation and visualize it according to different parameters.

These are the available functions:

- visualize the presence of RFI per frequency channel and baseline length.
- visualize the percentage flagged visibilities due to RFI per frequency channel. 
- visualize the increase in noice due to RFI per frequency channel.
- estimate the shape of the PSF after flagging of RFI.

check out the [WiKi](https://github.com/Fil8/RFInder/wiki) for a complete illustration of `RFInder`.

***
### Usage

All `rfinder` can be run in an automated way as `python rfi_pipeline.py <path_to_parameter_file.yml>`, or via an `IPython`
[notebook](https://github.com/Fil8/RFInder/blob/master/tutorials/T2_rfinder_automated.ipynb). Commands and options of `RFInder` must be specified in a `.yml` [parameter file](https://github.com/Fil8/RFInder/wiki/Parameter-file)

These [tutorials](https://github.com/Fil8/RFInder/tree/master/tutorials) show the different capabilities of `rfinder`.

***

### Installation

**Requisites**
- RFInder makes use of the most common `python` packages (e.g. `numpy`, `scipy`, `astropy`). 
- The parameter file is in `yaml` format, hence [`pyaml`](https://anaconda.org/anaconda/pyyaml), and [`json`](https://anaconda.org/conda-forge/json-c) packages should be installed,
- The `logging` module is used to print out warnings.
- `.gif` file of multiple plots can be created if `ffmpeg` is installed.
- `casacore` is utilized to open casa tables.
- **beam_shape** uses [`wsclean`] option `wcclean --psf-only`. Instructions to download and install `wsclean` [(Offringa et al. 2014)](https://arxiv.org/abs/1407.1943) can be found [here](https://sourceforge.net/projects/wsclean/).

**Insallation instructions**
To install from source clone this repository. From terminal type:

```
git clone https://github.com/Fil8/RFInder.git
```

Then run:

```
cd RFInder && pip install .
```

This package will soon be available on PYPI, allowing:

```
pip install rfinder
```

**License**

This project is licensed under the GNU General Public License v3.0 - see [license](https://github.com/Fil8/RFInder/blob/master/LICENSE.md) for details.


 ***
 <p>&copy <sub> Filippo M. Maccagni 2018 </sub></p>
