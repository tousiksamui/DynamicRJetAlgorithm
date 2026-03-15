#####################################
### Installtion Guide ###############
#####################################

### Prerequisite
- FastJet >= 3.1.0
  (Please see [FastJet](https://fastjet.fr/) for FastJet installation)

### Alternative 1.
Please use this alternative if the path for FastJet is properly set. Check whether `fastjet-config --help` command is executed properly.
  ```bash
  $ ./configure --fastjet-config=fastjet-config
  $ make
  $ make install
  ```

### Alternative 2.
Please use this alternative if write permission is not set but the user has sudo access.
  ```bash
  $ ./configure --fastjet-config=fastjet-config
  $ make
  $ sudo make install
  ```

### Alternative 3.
Please use this alternative if the path for FastJet is not set, but th einstallation directory of FastJet is known (e.g., /path/to/fastjet/bin/). Check whether `/path/to/fastjet/bin/fastjet-config --help` command is executed properly.
  ```bash
  $ ./configure --fastjet-config=/path/to/fastjet/bin/fastjet-config
  $ make
  $ make install
  ```

### After Installation
- DynamicRJetPlugin.hh and libDynamicRPlugin.so will be installed within the FastJet installation path.
- The usage is the same as the standard FastJet package.
- To use this plugin, please use the `-lDynamicRPlugin` flag.
- Check the examples directory for more examples.

### Citation
If you use the **Dynamic Radius Jet Clustering Algorithm** implemented in this repository, please cite the original algorithm paper:
- **B. Mukhopadhyaya, T. Samui, R. K. Singh**  
  *Dynamic Radius Jet Clustering Algorithm*  
  Journal of High Energy Physics **04** (2023) 019 [10.1007/JHEP04(2023)019](https://doi.org/10.1007/JHEP04(2023)019)

The algorithm has been used in several phenomenological studies in different collider environments. If relevant to your work, please consider citing the following papers:

**Studies with realistic LHC pile-up**
- **A. Ghosh, S. Ghosh, S. Mitra, T. Samui, R. K. Singh**  
  *Improving Sensitivity of Vectorlike Top Partner Searches with Jet Substructure*  
  Physical Review D **113** (2026) 035034 [10.1103/dhs6-gcn9](https://doi.org/10.1103/dhs6-gcn9)

**Zero pile-up studies**
- **B. Mukhopadhyaya, S. Samanta, T. Samui, R. K. Singh**  
  *Improving probes of hAA coupling in the Type-X two Higgs doublet model scenario: The crucial role of tau-jet charge identification*  
  Physics Letters B **856** (2024) 138949 [10.1016/j.physletb.2024.138949](https://doi.org/10.1016/j.physletb.2024.138949)

- **A. Ghosh, P. Konar, T. Samui, R. K. Singh**  
  *Jet substructure probe on scalar leptoquark models via top polarization*  
  Journal of High Energy Physics **07** (2025) 145 [10.1007/JHEP07(2025)145](https://doi.org/10.1007/JHEP07(2025)145)

### Note
- This module of dynamic radius jet algorithm is developed within the framework of the FastJet package.
- For more details on FastJet, please visit [FastJet](https://fastjet.fr/).

