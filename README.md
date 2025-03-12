# South Galapagos-1

## Code supporting the manuscript "Industrial Gamma Radiation Signals: An untapped resource to reconstruct Neogene paleoclimate off Northwest Australia"

Natural Gamma Radiation data used for this paper is stored in .csv files and is freely available through the National Offshore Petroleum Information Management System (NOPIMS, Geoscience Australia) (https://www.ga.gov.au/nopims). The core-based NGR data for the IODP Sites  is available through the IODP repository (https://web.iodp.tamu.edu/LORE/?appl=LORE&reportName=ngr). The scientific downhole logging data is available at https://mlp.ldeo.columbia.edu/logdb/. 

The repository contains the R scripts for the DTW calculations performed and presented for the Paleoceanography and Paleoclimatology paper. The separate R script ´Custom Step Pattern.R´ has the custom step pattern asymmetricP1.1 required for the DTW correlation technique, which must be run first. The `U1482-SG1_Seismic-NGR.R´ contains the DTW correlation seen in Figure S1. 

`SG1_Seismics-NF_Cyclostratigraphy.R´ file contains the actual tuning information, spectral and wavelet analyses, and the script to plot all the figures. 

`Picard1-Minilya1.R´ and `Picard1-U1464.R´ scripts for the DTW correlation must be run first before the `NGR-IndustrialSites.R´. This file contains the lines for plotting Figure 12. `U1447-U1501-U1482-SG1-NGR.R´ must be used for plotting Figure 13.

The updated CENOGRID data and metadata between 15.55 and 20 Ma is updated by replacing data from ODP Site 1264 (Liebrand et al., 2016) and IODP Site U1490 (Holbourn et al., 2020). The updated data and metadata is provided in the data folder.

## License

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

## Dependencies

- Giorgino T (2009). “Computing and Visualizing Dynamic Time Warping Alignments in
  R: The dtw Package.” _Journal of Statistical Software_, *31*(7), 1-24.
  doi:10.18637/jss.v031.i07 <https://doi.org/10.18637/jss.v031.i07>.
  
- Meyers, S.R. (2014). Astrochron: An R Package for Astrochronology.
  https://cran.r-project.org/package=astrochron

- Signorell A (2024). _DescTools: Tools for Descriptive Statistics_. R package
  version 0.99.54, <https://CRAN.R-project.org/package=DescTools>.

- Gouhier, T., Grinsted, A., & Simko, V. (2018). R package biwavelet: Conduct univariate and bivariate wavelet analyses. Version 0.20, 19. https://gitplanet.com/project/biwaveletGouhier, T., Grinsted, A., & Simko, V. (2018). R package biwavelet: Conduct univariate and bivariate wavelet analyses. Version 0.20, 19. <https://gitplanet.com/project/biwavelet>

- Callahan, J., Casey, R., Sharer, G., Templeton, M., & Trabant, C. (2022). IRISSeismic: Classes and Methods for Seismic Data Analysis (Version 1.6.6) [Computer software]. <https://cran.r-project.org/web/packages/IRISSeismic/index.html>

- Liebrand, D., Beddow, H. M., Lourens, L. J., Pälike, H., Raffi, I., Bohaty, S. M., Hilgen, F. J., Saes, M. J. M., Wilson, P. A., van Dijk, A. E., Hodell, D. A., Kroon, D., Huck, C. E., & Batenburg, S. J. (2016). Cyclostratigraphy and eccentricity tuning of the early Oligocene through early Miocene (30.1–17.1 Ma): Cibicides mundulus stable oxygen and carbon isotope records from Walvis Ridge Site 1264. Earth and Planetary Science Letters, 450, 392–405. <https://doi.org/10.1016/j.epsl.2016.06.007>

- Holbourn, A., Kuhnt, W., Kulhanek, D. K., Mountain, G., Rosenthal, Y., Sagawa, T., Lübbers, J., & Andersen, N. (2024). Re-organization of Pacific overturning circulation across the Miocene Climate Optimum. Nature Communications, 15(1), 8135. <https://doi.org/10.1038/s41467-024-52516-x>
