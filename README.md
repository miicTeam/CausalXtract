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
Concerning the immune cells, the parameters are:\
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

As the dynamic of your cells will likely differ from the ones used in the CausalXtract publication, there are two important parameters $\tau$ and $\delta\tau$ to tune for the causal discovery part on your data with tMIIC.
$\tau$ is the maximum number of layers back in time used to run the causal
discovery and $\delta\tau$ (1 by default) represents the number of timesteps between each layer.\
\
A possible way is to determine the $\tau$ parameter as follow:
- determine $L$, the minimum length of all the trajectories
- compute autocorrelation $acf_{t,v}$ for each trajectory $t$ and variable $v$ 
  over a maximum of $L$ time steps back in time\
- look for each trajectory and variable when the auto-correlation vanishes, retain half of the vanishing time to pick up the auto-correlation and compute the $\alpha_{t,v}$\
\
$l_{t,v} = \lfloor \frac{min\ (l\ |\ acf_{t,v}[l]\ <\ 0.05)}{2} \rceil$ \
$\alpha_{t,v} = acf_{t,v}[l_{t,v}]^{\frac{1}{l_{t,v}}}$ \
\
The $\alpha_{t,v}$ are averaged over the $T$ trajectories and $V$ variables and the mean is used to determine the $\tau$ parameter as twice the relaxation time in order to be sure to cover the complete dynamic of the system.
\
$\overline{\alpha} = \frac{\sum_{t=1}^{T} \sum_{v=1}^{V} \alpha_{t,v}}{T\ *\ V}$ \
$\tau = \lfloor 2 * \frac{1\ +\ \overline{\alpha}}{1\ -\ \overline{\alpha}}\rceil$ \
If the $\tau$ parameter leads to too much nodes in the time lagged graph, 
you can increase the $\delta\tau$ parameter and reduce the $\tau$ so that $\tau_{initial} \simeq \tau_{final} * \delta\tau_{final}$. 
As example, in the CausalXtract publication, the parameters used were 
$\tau=11$ and $\delta\tau=5$ resulting in 182 nodes in the inferred graph, still 
achievable on a recent computer.

You can find more information one the causal discovery module by calling the documentation of the main function `?miic` from R.

Example of tMIIC on a toy dataset:
```R
library(miic)

# EXAMPLE: TOY MODEL COVIDCASES
data(covidCases)
# execute tMIIC
# supplying the tau parameter with a value >= 1 switches miic into temporal mode 
# here, we perform a temporal causal discovery with 2 time steps back in history 
tmiic.res <- miic(input_data = covidCases, tau = 2)

# to plot the default graph (compact)
if(require(igraph)) {
 plot(tmiic.res)
}
# to plot the full temporal network
if(require(igraph)) {
  plot(tmiic.res, display="lagged")
}
```

## Documentation

You can find the documentation pages in the "help" folder of the tMIIC module, or use the R functions `help()` and `?`.

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
