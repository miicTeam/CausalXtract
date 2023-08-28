# CausalXtract

This repository contains the source code for CausalXtract, a flexible pipeline to extract causal effects from live-cell time-lapse imaging data.

CausalXtract is composed of two parts: the cell feature module **CellHunter+** and the causal discovery module **tMIIC**:

* CellHunter+ is a method that segments cells, tracks them and extracts features. The segmentation is based on CHT (Circular Hough Transform); the tracking is based on the Munkres algorithm. 
It extends the functionalities of CellHunter (that implements segmentation and tracking) by adding the feature extraction module.

* tMIIC is the temporal version of MIIC (**M**ultivariate **I**nformation based **I**nductive **C**ausation), a method based on constraint-based approaches that learns a large class of causal or non-causal graphical models from purely observational data while including the effects of unobserved latent variables.


## References

CausalXtract, submission pending.

Nguyen M., De Ninno A., Mencattini A., Mermet-Meillon F., Fornabaio G., Evans SS., Cossutta M., Khira Y., Han W., Sirven P., Pelon F., Di Giuseppe D., Bertani FR., Gerardino A., Yamada A., Descroix S., Soumelis V., Mechta-Grigoriou F., Zalcman G., Camonis J., Martinelli E., Businaro L., Parrini MC.; Dissecting Effects of Anti-cancer Drugs and Cancer-Associated Fibroblasts by On-Chip Reconstitution of Immunocompetent Tumor Microenvironments; Cell Rep.; 2018 Dec 26;25(13):3884-3893.e3. [doi: 10.1016/j.celrep.2018.12.015](https://doi.org/10.1016/j.celrep.2018.12.015)

Cabeli V., Verny L., Sella N., Uguzzoni G., Verny M., Isambert H.; Learning clinical networks from medical records based on information estimates in mixed-type data; PLoS computational biology., 2020. [doi:10.1371/journal.pcbi.1007866](https://doi.org/10.1371/journal.pcbi.1007866) | [code](https://github.com/vcabeli/miic_PLoS)

Verny L., Sella N., Affeldt S., Singh PP., Isambert H.; Learning causal networks with latent variables from multivariate information in genomic data;  PLoS Comput. Biol., 2017. [doi:10.1371/journal.pcbi.1005662](https://doi.org/10.1371/journal.pcbi.1005662)


## Prerequisites

CellHunter+ was written in **MATLAB R2022b**. Additional toolboxes required: **"Image Processing Toolbox"**, **"Statistics and Machine Learning Toolbox"**. The former is needed for both the "segmentation and tracking" module and the "feature extraction" module. The latter is needed only for the "feature extraction" module. CellHunter+ should work with other versions of MATLAB as well, but be aware that the definition of the features extracted could slightly change. Therefore, the resulting data obtained with CellHunter+ could be slightly different from the data obtained for the article.

As the tMIIC part contains R and C++ sources, you will need **R** and a compiler with support for **C++14** language features. tMIIC imports the following R packages: **ppcor, scales, stats, Rcpp** and, if you use the plotting feature, **igraph**. The R version used to develop tMIIC is the 4.0.5, so any R version equal or above 4.0.5 should be usable.

## Installation

To replicate the CausalXtract repository on your computer, you can use the **git clone** command:

```bash
git clone https://github.com/miicTeam/CausalXtract.git
```

From the R environment, you can install the causal discovery module tMIIC **remotely from github** or (after cloning the CausalXtract repository) **as an usual R source package**.
```R
# Remotely from github:
remotes::install_github ("miicTeam/CausalXtract/tMIIC")
# Locally from a clone of CaulsaXtract, i.e. from the CausalXtract folder:
install.packages ("./tMIIC", repos=NULL, type="source") 
```

## Quick start

A demo is available in the Demo folder: **CausalXtract_Demo.Rmd**.\
This R Markdown Notebook file allows you to run the entire pipeline described in the CausalXtract publication.

## Going further

To go further with CellHunter+, "the main_detection.m" file implements the "segmentation and tracking" module. Each ROI (Region of Interest) is a cropped video, obtained from the original video (reference to dataset: https://doi.org/10.5281/zenodo.7755700). The MCC (Main Cancer Cell) is placed at the centre of the crop and it is possible to observe CAFs, immune cells and other cancer cells. 
It is possible to change the file path to the specific videos path that the user wants to analyse.\
The video of the ROI, saved as .mat file, is a matrix where the third dimension represents the time, i.e. the number of frames. 
The outputs are the trajectories of the MCC and the trajectories of the immune cells.
An additional step, implemented in the "main_division_detection.m" file, allows to correct the "flickering" of the MCC's trajectory when it divides.
Finally, the "main_features.m" file implements the "feature extraction" module, in which the trajectories of the MCC and those of the immune cells are used to compute the features of interest for each ROI.
\
As the dimensions of your cells will likely be different from the ones used in the CausalXtract publication, you may need to modify the parameters in the "parameters_CellHunterPlus.csv" file.
1. The first parameter, "polarity", must be set to "bright" if bright cells are identified in a dark background. Otherwise, it must be set to "dark".\
Concerning the immune cells, the parameters are:
2. r_sp, the theoretical radius for detecting immune cells; 
3. Rmax_sp, the maximum distance for tracking immune cells, i.e. for linking two presumed instances of the same immune cell in two different frames to construct the trajectory of that immune cell; 
4. DP_sp, the number of frames after which the trajectory of an immune cell is stopped if the immune cell is not detected for that specific number of frames;
5. L_sp, the minimum length of immune cells trajectories that are returned as output in the tracking refining process;
6. r_std, a parameter that allows to delete the trajectories of presumed immune cells that do not move enough to be considered as such. If r_std is increased, less immune cells are considered.\
Concerning the cancer cells, the parameters are: 
7. r_tu; 8. Rmax_tu, 9. DP_tu, which are the analagous parameters of those for the immune cells defined above.
Additionally, 10. dist_tu is a parameter that imposes to detect only the cancer cells whose centre is less than dist_tu pixels away from the centre of the ROI.\
As example, in the CausalXtract publication, the parameters used are the following:
polarity="bright";
r_sp=4; Rmax_sp=20; D_sp=10; L_sp=10; r_std=4
r_tu=14; Rmax_tu=40; D_tu=70; dist_tu=30

Even if the dynamic of your cells will likely differ from the one in the CausalXtract publication, the causal discovery part tMIIC includes an automatic estimation of the temporal dynamic and will adapt accordingly.

To perform a first try on the features extracted by CellHunter+, the minimal parameters of tMIIC are:
```R
miic_res <- miic(input_data=dataframe_features, mode="TS")
```
It will produce a final (lagged) graph with only 50 nodes by default in order to run quickly, so you can have a first look on the result of tMIIC:
```R
# default compact plot
plot(miic_res)
# lagged plot
plot(miic_res, display="lagged")
```

To go even further, tMIIC has several parameters useful to know:
- max_nodes: can be used to increase (or decrease) the maximum number of nodes
in the lagged graph. tMIIC uses this parameter, once the temporal dynamic has been estimated, to compute the number of layers and the number of time steps between the layers in a way that covers the dynamic. The more nodes you allow in the final graph, the more time tMIIC will need to perform the discovery, but more precise or complete will be the discovery: more layers and/or smaller number of time steps between each layer. On recent computers or servers, values up to 200 or 300 are possible (depending on the number of time steps in your dataset). 
n_layers and delta_t: you can specify your own parameters for the number of layers and the number of time steps between layers. In such case, tMIIC will not perform the automatic dynamic estimation and will use the provided parameters instead. 
- state_order: the state_order is an optional data frame that can be used to specify extra characteristics per variable. The information that you are the more likely to be interested to supply are var_names, levels_increasing_order and is_contextual.
  * var_names: mandatory, is the name of the variables in the input dataset
  * levels_increasing_order: optional, can be used to specify an order for the discrete variables. A typical example is for logical variables as "Treatment", where we can add 0,1 as levels_increasing_order to display colored edges highlighting the negative and positive correlations with "Treatment".
  * is_contextual: optional, is frequently used for experimental conditions, that are set up from the start of each experiment and don't change over time. Values can be 0 or 1, where 0 is a normal variable and 1 indicates a contextual one. Variables defined as contextual can not be the consequence of any other variables in the dataset.

To have an example on how to set up a state_order, you can have a look on the one used in the CausalXtract publicaton that is in the demo folder: **CausalXtract_Publication_State_Order.tsv**. 
  
More information about the tMIIC parameters is also available by calling the documentation of the miic R package.

## Documentation

You can find the documentation pages about the tMIIC module by using the R functions `help()` and `?`.

'''R
help(miic)
?miic()
'''

## Authors
- Franck Simon
- Maria Colomba Comes
- Tiziana Tocci
- Louise Dupuis
- Vincent Cabeli
- Nikita Lagrange 
- Arianna Mencattini
- Maria Carla Parrini
- Eugenio Martinelli
- Hervé Isambert

## License
GPL-3
