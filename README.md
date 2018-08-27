# RFInder


This is a set of tools that have been developed in preparation of the [Apertif surveys](
https://www.astron.nl/astronomy-group/apertif/science-projects/apertif-science-projects).

The main function of `rfinder` is to identify the presence of RFI in an observation and visualize it according to different parameters.

These are the available functions:

- visualize the presence of RFI per frequency channel and baseline lenght.
- visualize the percentage flagged visibilities due to RFI per frequency channel. 
- visualize the increase in noice due to RFI per frequency channel.
- estimate the shape of the PSF after flagging of RFI.

`rfinder` is run using a `.yml` [parameter file](https://github.com/Fil8/RFInder/wiki/Parameter-file) as `python rfinderpipeline.py <path_to_parameter_file.yml>`, or through a `IPython`
[notebook](https://github.com/Fil8/RFInder/blob/master/tutorials/T2_rfinder_automated.ipynb). 

These [tutorials](https://github.com/Fil8/RFInder/tree/master/tutorials) can guide you through the different capabilities of `rfinder`.

***

### Installation

**Requisites**
- RFInder makes use of the most common `python` packages (e.g. `numpy`, `scipy`, `astropy`). 
- The parameter file is in `yaml` format, hence [`pyaml`](https://anaconda.org/anaconda/pyyaml), and [`json`](https://anaconda.org/conda-forge/json-c) packages should be installed,
- Tutorials make use of [`tabulate`](https://pypi.org/project/tabulate/) and [`glob`](https://anaconda.org/conda-forge/glob2) for fancy outputs.
- **beam_shape** uses [`wsclean`] option `wcclean --psf-only`. Instructions to download and install `wsclean` [(Offringa et al. 2014)](https://arxiv.org/abs/1407.1943) can be found [here](https://sourceforge.net/projects/wsclean/).

**Insallation instructions**
- Clone this repository. From terminal type:

```
git clone https://github.com/Fil8/RFInder.git
```

- add `rfinder` directory to `PYTHONPATH` (for permanent use save it in your `.cshrc_profile`, or `.bash_profile`, respectively)

```
setenv PYTHONPATH $path_to_rfindser:${PYTHONPATH}

export PYTHONPATH=$PYTHONPATH:path_to_rfinder
```

- change path at `line16` of `rfinder.py` as follows: `sys.path.append('/path-to-rfinder/RFInder_modules/')` 
 
 ***
 <p>&copy <sub> Filippo M. Maccagni 2018 </sub></p>
