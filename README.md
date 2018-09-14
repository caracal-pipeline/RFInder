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

To call `rfinder` automatically with the parameters stored in `rfinder_default.yml` add the following lines to your `.cshrc/.bashrc` file:

```
setenv RFI '<path_to_rfinder_pipeline.py>'
alias rfinder 'python $RFI/rfi_pipeline.py $RFI/rfinder_default.yml'
```

These [tutorials](https://github.com/Fil8/RFInder/tree/master/tutorials) can guide you through the different capabilities of `rfinder`.

***

### Installation

**Requisites**
- RFInder makes use of the most common `python` packages (e.g. `numpy`, `scipy`, `astropy`). 
- The parameter file is in `yaml` format, hence [`pyaml`](https://anaconda.org/anaconda/pyyaml), and [`json`](https://anaconda.org/conda-forge/json-c) packages should be installed,
- The `logging` module is used to print out warnings.
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

### Functions


We developed a small software that reading an MS file can show the Radio Frequency Interference (RFI) in different ways, RFInder (LINK). This tool can be useful to assess the quality of a radio interferometric observation and/or the amount of flagged visibilities of an observation. 

The main function of RFInder is to sort the visibilities (or the flags) in an MS file by baseline, frequency and time, with the baselines sorted by length. This allows show the amount of flagged visibilities of an observation for each baseline and frequency channel and highlight the dependence of RFI by the length of the baseline receiving that signal. RFInder in  to investigates the RFI signal of different correlations (e.g. XX,YY,XY,YX) in two main ways. One before running an automated flagger (RFI investigation) and one for an a posteriori investigation of the flagged visibilities. 

Two main options are available:
    
- RFI investigation: for each baseline the r.m.s is measured in a RFI free frequency interval selected by the user. For each baseline, is considered as RFI all signal that has amplitude > (threshold x r.m.s). The threshold is selected by the user (= 5 in the Figures below). For each baseline, the percentage of visibilities that are RFI is recorded. This, considered as RFI signal, is plotted with the baselines ordered by length, as shown in the Figure below.

- FLAGS investigation: RFInder opens the FLAG column of an MS file (the flags made by the observer, and/or the result of an automated flagger) and plots the percentage of flags for each baseline. Baselines are ordered by length. This tool is useful to see how many visibilities have been flagged in an observation. The Figure below shows the FLAGS of a run of aoflagger on the same observation of the Figure above. 

RFInder can select all visibilities and visibilities above and below a length given by the user. RFInder averages the %of RFI (or Flags) of those baselines and estimates how the natural weighted r.m.s. of the observation given the loss of visibilities. Figure A shows the % of visibilities for baselines above and below 200 m at the MeerKAT telescope. 

RFInder can also divide an observation in small time sets of length given by the user. Selecting a spectral window (of width given by the user) we can look at the % of RFI by altitude and azimuth of an observation. Figure XXXX shows the % of RFI every 10 minutes an observation of XXX hours in the spectral window XXXX.


 ***
 <p>&copy <sub> Filippo M. Maccagni 2018 </sub></p>
