The present code is my implementation of the seismic data processing algorithm described by G. D. Bensen, M. H. Ritzwoller, M. P. Barmin, A. L. Levshin, F. Lin, M. P. Moschetti, N. Shapiro, Y. Yang in "Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements". It also includes several useful tools for visualizing the results of this processing. 

Initially I wrote this code to work with the dataset that ultimately became the basis of the paper "Crustal structure beneath Central Kamchatka inferred from ambient noise tomography" by I. Egorushkin, I. Koulakov, A. Jakovlev, H.-H. Huang, E. I. Gordeev, I. Abkadyrov, D. Chebrov. Now I aim to rework it into a comprehensive package that can be easily customized for any dataset.

A couple of remarks about the current structure of the package:
- the core functions are stored in the folder with the same name;
- the main processing steps are performed in scripts located in the root directory;
- 'other' contains a number of situationally useful scripts. 
