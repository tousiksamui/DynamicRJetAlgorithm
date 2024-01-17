#####################################
### Installtion Guide ###############
#####################################

### Prerequisite
- FastJet >= 3.3.2
  (Please see [FastJet](https://fastjet.fr/) for FastJet installation)

### Alternative 1.
Please use this alternative if path for FastJet is properly set. Check whether `fastjet-config --help` command is executed properly.
  $ ./configure --fastjet-config=fastjet-config
  $ make
  $ make install

### Alternative 2.
Please use this alternative if write permission is not set but the user has sudo access.
  $ ./configure --fastjet-config=fastjet-config
  $ make
  $ sudo make install

### Alternative 3.
Please use this alternative if path for FastJet is not set but installation directory of FastJet is known (e.g., /path/to/fastjet/bin/). Check whether `/path/to/fastjet/bin/fastjet-config --help` command is executed properly.
  $ ./configure --fastjet-config=/path/to/fastjet/bin/fastjet-config
  $ make
  $ make install


### After Installation
- DynamicRJetPlugin.hh and libDynamicRPlugin.so will be installed within the FastJet installation path.
- The usage is the same as the standard FastJet package.
- To use this plugin, please use the `-lDynamicRPlugin` flag.
- Check the example directory for more examples.


### Note
- This module of dynamic radius jet algorithm is developed within the framework of the FastJet package.
- For more details on FastJet, please visit [FastJet](https://fastjet.fr/).
- The implementation is based on [DOI: 10.1007/JHEP04(2023)019](https://doi.org/10.1007/JHEP04(2023)019).
- If you use this for your work, please cite this reference and the reference to FastJet package.

