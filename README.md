# Eddy Tracking South Atlantic Bight

This folder contains the codes and outputs for the "Proceedings of The Royal Society A" paper of D Lindo-Atichati, K Abbot, J Gula , and JC McWilliams (manuscript in preparation)  

1) We ran the last version of py-eddy-tracker on Bitbucket (https://bitbucket.org/emason/py-eddy-tracker, update 11/30/2018). That is a major change made by the developer Evan Mason in 2018. Specifically, we used the module dataset.UnRegularGridDataset. This new module and version finally enabled us to accurately process the rotated, irregular grid dataset of the ROMS ocean model downscaled for the North Western Atlantic (Gula et al. 2016). That code produced eddy identification and properties in the region for 18 years. The outputs are in the files Anticyclonic.nc and cyclonic.nc.

2) We binned oceanographic outputs (e.g. vorticity), eddies, and eddy properties into a 1 deg x 1 deg grid for the North Atlantic Bight. Binning codes were written by Katy Abbot, the marvelous Helen Fellow I mentor at the American Museum of Natural History.
