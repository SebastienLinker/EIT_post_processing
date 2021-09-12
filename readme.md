# Post processing Electrical Impedance Tomography

Electrical impedance tomography (EIT) with neural networks.

This repo contains the source code to reproduce the results from those papers

> [**A Post-Processing Method for Three-Dimensional Electrical Impedance Tomography**](https://doi.org/10.1038/s41598-017-07727-2),
> Martin, S., Choi, C.T.M. A Post-Processing Method for Three-Dimensional Electrical Impedance Tomography. Sci Rep 7, 7212 (2017) 
> *https://doi.org/10.1038/s41598-017-07727-2*

> [**A novel post-processing scheme for two-dimensional electrical impedance tomography based on artificial neural networks*](ttps://doi.org/10.1371/journal.pone.0188993),
> Martin, S., Choi, C.T.M. A novel post-processing scheme for two-dimensional electrical impedance tomography based on artificial neural networks. PLOS ONE 12(12): e0188993
> *ttps://doi.org/10.1371/journal.pone.0188993*

## How to use

Download EIDORS from the original SVN repo: http://eidors3d.sourceforge.net/

Download MatGeom from Github (included as a submodule): https://github.com/mattools/matGeom

Install Netgen: https://ngsolve.org/downloads

Add EIDORS and MatGeom into your Matlab path, add this repository into your Matlab path as well

The `demo` directory contains some sample codes to reproduce the results. You can update the parameters at the beginning of the file

### 2D EIT

You can choose whether to add noise, distortion, the target type, them the script train the models, generate some additional random targets and displays the results

### 3D EIT

Similar to 2D EIT. Note that you may need a huge amount of memory to solve some large problems. Finite elements models with a smaller number of nodes have been designed for rapid testing, but those models won't give any good result due to coarsity of the models.

## Requirements

### MATLAB 

This code was originally developed for different versions of Matlab, from R2013a from R2016a, and successfully tested again under version R2020b

### MATLAB Toolboxes

Required:
Neural networks toolbox

Parallel computing toolbox can make things faster, although it was crashing on latest Matlab version
Communications toolbox and Symbolic math toolbox are required to use the noise functions

### EIDORS

EIDORS versions 3.8 to 3.10 are working correctly

### MatGeom

Version 1.2.2

### Netgen

Tested with Netgen 6.2 https://ngsolve.org/downloads

