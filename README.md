# RFInder

**Insallation instructions**

```
pip install rfinder
```

To create a local repository, type:

```
git clone git@github.com:caracal-pipeline/RFInder.git
```

***

**Requisites**

For a successfull installation make sure to have installed the following packages.

- RFInder makes use of the most common `python` packages (e.g. `numpy`, `scipy`, `astropy`). 
- The parameter file is in `yaml` format, hence [`pyaml`](https://anaconda.org/anaconda/pyyaml), and [`json`](https://anaconda.org/conda-forge/json-c) packages should be installed,
- The `logging` module is used to print out warnings.
- `.gif` file of multiple plots can be created if `ffmpeg` is installed.
- `casacore` is utilized to open casa tables.
    - install it with `python_casacore`: `pip install python-casacore` or `conda install -c conda-forge python-casacore` 
- `texmaker` to plot latex fancy formulae
- `dvipng`
- `python tk`
- RFInder class `rfinder_uzero.py` which can be used to solve the [u=0 problem](https://archive-gw-1.kat.ac.za/public/repository/10.48479/bhpj-nz95/index.html) makes use of `wsclean`, which must have been previously installed.


***
**Description**

This is a set of tools that have been developed in preparation of the Apertif & MeerKAT surveys.

The main function of `rfinder` is to identify the presence of RFI in an observation and visualize it according to different parameters. Two are the main functions:

- estimate the RFI present in an MS file through a sigma clipping (`rms_clip`)
- read the `FLAG` column of an MS file (`use_flags`) and summarize how RFI affects the data products of an observation. 

These are the products that `rfinder` provides and summarizes in an `.html` file:

- presence of RFI per frequency channel and baseline length.
- percentage flagged visibilities due to RFI per frequency channel. 
- increase in noice due to RFI per frequency channel.
- estimated noise per frequency channel, assuming natural weighting. 

check out the [WiKi](https://github.com/Fil8/RFInder/wiki) for a complete illustration of `RFInder`.

***
**Usage**

RFInder takes its variables from a default parameter file and from terminal, if any are given. 

From your current working directory typying `rfinder` this message will be shown: 

```
------ Reading default installation parameter file ------

MSNAME & telescope missing
              		please edit rfinder_default.yml in your current directory
              		or run: rfinder -i msname -fl <num> -tel <meerkat,apertif,wsrt>
              		(assuming the observation is located in your current directory)
                    

------ RFInder out ------
```

Hence, you have to set the name of the MSfile you wish to analyse. There are two ways to do this. By specifying from terminal the path to the msfile from your current directory, the field number of the source you whish to analyse, and the telescope of the observation:

```
 rfinder -i msname -fl <num> -tel <meerkat,apertif,wsrt>
```

or, editing the `rfinder_default.yml` configuration file that has been copied in your current directory (workdir, in the configuration file). 

This configuration file is read automatically by RFInder through the command `rfinder`. A short explanation of the parameters is given in the configuration file, and by typing `rfinder -h` (see below).

If you wish to use a different configuration file (at your own risk!!), type: `rfinder -c <path_to_configuration_file>`.

**Minimal instructions**

- Default `rfinder` will scan the MSfile in chunks of 10 minutes averaging 10 channels together. The output product will be an `html` file where the `gis` scan through the time steps to show the identified RFI/flags.

- Running `rfinder -noCh` after `rfinder` will produce a `full_report.html` file containing both the analysis over time steps and the analysis of the dataset as a whole.

- Running `rfinder -noCh -noMov` will analyse the full dataset as a whole and generate the `full_report.html` without embedded movies.

_Attention_: the option `rfinder -noCh` will end with a report successfully generated, only if it is run after `rfinder`. Otherwise run `rfinder -noCh -noMov`.

(These [tutorials](https://github.com/Fil8/RFInder/tree/master/tutorials) show the different capabilities of `rfinder`. **outdated**)


**Output products**

If `rfinder` runs correctly, you will find the following output products in your current directory: 

- the folder `rfi_pol` in your current directory, or in the directory specified by the `-odir` parameter (`pol` is the stokes parameters for which you analysed RFI). 
	- Within, there are the `.html` reports that you wished to generate. 
- The configuration file `rfinder_default.yml` contains the parameters of the last run.
- A `log` of the commands run by the program is stored in `log-rfinder.log`, in your working directory.

**Help**

`rfinder -h` will show you a (minimal) help:

```
RFInder: package to visualize the flagged RFI in a dataset
version 1.2.0
install path /path/to/RFInder/rfinder
Filippo Maccagni <filippo.maccagni@gmial.com>

Usage:  [OPTIONS]

Options:
  -idir, --input-dir Directory    Full path to the working directory
  -c, --config File               RFInder configuration file (YAML format)
  -i, --msname MS                 Name of the input Measurement Set (MS) file
  -fl, --field int                Field ID of the target in the file
  --cleanup-enable / --no-cleanup-enable
                                  Remove intermediate results
  -j, --ncpu int                  Number of CPUs to get the total flagged data
  -l, --label str                 Label of the output directory. Result:
                                  rfi_<stokes>_<label>
  -tel, --telescope-name str      Name of the telescope.
  --telescope-diameter float      Diameter of the telescope in meters
  --telescope-tsyseff float       Effective system temperature in Kelvin
  --telescope-long float          Longitude of the telescope
  --telescope-lat float           Latitude of the telescope
  --telescope-height float        Height/altitude of the telescope location
  -rfi, --rfi-enable / --no-rfi-enable
                                  Enable or disable the RFI detection module
  -pol, --polarization str        Polarization type (e.g., 'xx', 'yy', 'q')
  --bad-antenna str,str,...       List of bad antennas
  -mode, --rfimode str            Mode of RFI detection ('rms_clip' or
                                  'use_flags')
  --sigma-clip, --sig float       Threshold for RFI identification
  -fint, --frequency-interval float,float,...
                                  Frequency range to measure average STD of
                                  visibilities
  --baseline-cut int              Cutoff baseline length
  --chunks-time-enable / --no-chunks-time-enable
                                  Enable splitting by time intervals
  -tStep, --chunks-time-step int  Time chunk size in minutes
  --chunks-spw-enable / --no-chunks-spw-enable
                                  Enable splitting by spectral windows
  -spwAv, --chunks-spw-width int  Channel width of rebinned output table in
                                  MHz
  --plot-details-enable / --no-plot-details-enable
                                  Enable detailed plotting
  --plot-details-plot-noise str   Type of noise/RFI to plot ('rfi', 'noise',
                                  or 'noise_factor')
  --plot-details-plot-long-short / --no-plot-details-plot-long-short
                                  Plot all baselines or only long/short
                                  baselines
  --plot-details-plot-eps / --no-plot-details-plot-eps
                                  Generate EPS plots
  --plot-details-movies-2d-gif / --no-plot-details-movies-2d-gif
                                  Generate 2D GIF movies
  --plot-details-movies-1d-gif / --no-plot-details-movies-1d-gif
                                  Generate 1D GIF movies
  --plot-details-movies-altaz-gif / --no-plot-details-movies-altaz-gif
                                  Generate Alt/Az GIF movies
  --plot-details-movies-in-report / --no-plot-details-movies-in-report
                                  Include movies in the generated report
  --plot-summary-enable / --no-plot-summary-enable
                                  Enable summary plotting of % flagged
                                  visibilities
  --plot-summary-axis str,str,...
                                  Axes for summary plotting
  --plot-summary-antenna str      Select an antenna for summary
  --plot-summary-freq-bin int     Bin frequency channels for summary plotting
  --plot-summary-report / --no-plot-summary-report
                                  Generate HTML report with results
  -odir, --output-dir Directory   Full path to the working output directory
  -h, --help                      Show this message and exit.
  --version                       Show the version and exit.

Run a command. This can be:

rfinder 
rfinder -c path_to_config_file.yml
rfinder -i <ngc1399.ms> -fl <num> -tel <meerkat/apertif/wsrt>
rfinder -i <ngc1399.ms> -fl <num> -tel <meerkat/apertif/wsrt> -rfi -mode rms_clip
```

***

**License**

This project is licensed under the GNU General Public License v3.0 - see [license](https://github.com/Fil8/RFInder/blob/master/LICENSE.md) for details.


 ***
 <p>&copy <sub> Filippo M. Maccagni 2018 </sub></p>
