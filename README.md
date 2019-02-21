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

The main function of `rfinder` is to identify the presence of RFI in an observation and visualize it according to different parameters.

These are the available functions:

- visualize the presence of RFI per frequency channel and baseline length.
- visualize the percentage flagged visibilities due to RFI per frequency channel. 
- visualize the increase in noice due to RFI per frequency channel.
- estimate the shape of the PSF after flagging of RFI.

check out the [WiKi](https://github.com/Fil8/RFInder/wiki) for a complete illustration of `RFInder`.

***
**Usage**

`rfinder rfinder_default.yml`

Commands and options of `RFInder` are specified in the yaml file

These [tutorials](https://github.com/Fil8/RFInder/tree/master/tutorials) show the different capabilities of `rfinder`.

***

**License**

This project is licensed under the GNU General Public License v3.0 - see [license](https://github.com/Fil8/RFInder/blob/master/LICENSE.md) for details.


 ***
 <p>&copy <sub> Filippo M. Maccagni 2018 </sub></p>
