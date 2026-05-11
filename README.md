# morphoXS

This repository contains code to extract topography along a channel flowline (e.g. thalweg or centerline) as well as at channel-perpendicular cross-sections:

![Demo](https://github.com/shelbyahrendt/morphoXS/blob/main/demo_fig.png?raw=true)


## Installation

Install with conda:

1. Clone the repository to your local machine: ```git clone https://github.com/shelbyahrendt/morphoXS.git```
2. Within the repository tree run the following in your terminal ```conda env create -f environment.yml``` (this may take several minutes)
3. Activate the environment ```conda activate morphoXS```

## File Organization

* ```transect_extraction_functions.py``` contains functions to place equidistant points along a channel line, generate transects perpindicular to these points, as well as calculate channel curvature
* ``` demo_transect_extraction.ipynb``` demonstrates to how apply functions to example valley-centerline shapefile located in ```data``` as well as extract long-profile and cross-sectional topography

## Run Demo

1. Open Jupyter Lab (your morphoXS environment must be activated)
2. Select ```demo_transect_extraction.ipynb``` from the file tree to open

