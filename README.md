# RFInder

**Insallation instructions**

```
pip install rfinder
```

To create a local repository, type:

```
git clone https://github.com/Fil8/RFInder
```

***

**Requisites**

For a successfull installation make sure to have installed the following packages.

- RFInder makes use of the most common `python` packages (e.g. `numpy`, `scipy`, `astropy`). 
- The parameter file is in `yaml` format, hence [`pyaml`](https://anaconda.org/anaconda/pyyaml), and [`json`](https://anaconda.org/conda-forge/json-c) packages should be installed,
- The `logging` module is used to print out warnings.
- `.gif` file of multiple plots can be created if `ffmpeg` is installed.
- `casacore` is utilized to open casa tables.
- `texmaker` to plot latex fancy formulae
- `dvipng`
- `python tk`



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
usage: rfinder [-h] [-v] [-c CONFIG] [-w WORKING_DIR] [-odir OUTPUT_DIR]
               [-i INPUT] [-fl FIELD] [-tel TELESCOPE] [-mode RFIMODE]
               [-pol POLARIZATION]
               [-fint [FREQUENCY_INTERVAL [FREQUENCY_INTERVAL ...]]]
               [-spwAv SPW_AV] [-tStep TIME_STEP] [-sig SIGMA_CLIP]
               [-baseCut BASELINE_CUT] [-noCh] [-yesCh] [-noSpw] [-yesSpw]
               [-noClp] [-yesClp]

RFInder: package to visualize the flagged RFI in a dataset

version 1.0.3

install path /home/maccagni/programs/RFInder/rfinder

Filippo Maccagni <filippo.maccagni@gmial.com>

optional arguments:
  -h, --help            Print help message and exit
  -v, --version         show program's version number and exit
  -c CONFIG, --config CONFIG
                        RFInder configuration file (YAML format)
  -idir INPUT_DIR, --input_dir WORKING_DIR
                        select working directory (MS file assumed to be here)
  -odir OUTPUT_DIR, --output_dir OUTPUT_DIR
                        select output directory
  -i INPUT, --input INPUT
                        input ['MS'] file
  -fl FIELD, --field FIELD
                        select field of MS file to analyze
  -tel TELESCOPE, --telescope TELESCOPE
                        select telescope: meerkat, apertif, wsrt
  -mode RFIMODE, --rfimode RFIMODE
                        select mode where to investigate RFI: use_flags or
                        rms_clip
  -pol POLARIZATION, --polarization POLARIZATION
                        select stokes parameter: xx, yy, xy, yx, q (also in
                        CAPS)
  -fint [FREQUENCY_INTERVAL [FREQUENCY_INTERVAL ...]], --frequency_interval [FREQUENCY_INTERVAL [FREQUENCY_INTERVAL ...]]
                        select frequency interval where to measure noise in
                        GHz
  -spwAv SPW_AV, --spw_av SPW_AV
                        select number of channels to average
  -tStep TIME_STEP, --time_step TIME_STEP
                        select time step in minutes in which divide the
                        analysis of the MSfile
  -sig SIGMA_CLIP, --sigma_clip SIGMA_CLIP
                        select sigma clip for rms_clip mode to find RFI
  -baseCut BASELINE_CUT, --baseline_cut BASELINE_CUT
                        select cut in baseline lenght [m] for differential RFI
                        analysis
  -noCh, --no_chunks    desable chunking in time
  -yesCh, --yes_chunks  enable chunking in time
  -noSpw, --no_spw_av   desable averaging in channels
  -yesSpw, --yes_spw_av
                        enable averaging in channels
  -noClp, --no_cleanup  desable cleanup of intermediate products
  -yesClp, --yes_cleanup
                        enable cleanup of intermediate products

Run a command. This can be: 
rfinder 
rfinder -c path_to_config_file.yml
rfinder -i <ngc1399.ms> -fl <num> -tel <meerkat/apertif/wsrt>
```

***

**License**

This project is licensed under the GNU General Public License v3.0 - see [license](https://github.com/Fil8/RFInder/blob/master/LICENSE.md) for details.


 ***
 <p>&copy <sub> Filippo M. Maccagni 2018 </sub></p>
