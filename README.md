# BTRRwGGM
Bayesian Tensor Response Regression with a Gaussian Graphical Model

## Description
This is a small package that allows for the data generation and analysis following the manuscript 
by Spencer, Guhaniyogi, and Prado, which has been submitted to Biometrics for publication.

## Model
The Bayesian model formulation takes *G* different tensors of dimension *D*+2 as a 
response variable in a regression setting:

<img src="https://latex.codecogs.com/svg.latex?\mathbf{Y}_{g,t,i}&space;=&space;\mathbf{B}_{g}x_{t,i}&space;&plus;&space;d_{g,i}&space;&plus;&space;\mathbf{E}_{g,t,i}." title="\mathbf{Y}_{g,t,i} = \mathbf{B}_{g}x_{t,i} + d_{g,i} + \mathbf{E}_{g,t,i}." />
  
In the context of an fMRI scan, if one is examining a slice of the brain, then *G* would denote the number of mutually 
exclusive regions of interest, and *D* would be equal to 2. The two additional dimensions on the response tensor 
correspond to the time *t* of the scan, and the subject *i*, respectively. The variable *d*<sub>*g,i*</sub> is a 
subject-region effect, which should be centered around zero after preparing the fMRI data (details below). This model 
assumes that the elements in the error tensor, **E**<sub>*g,i,t*</sub>, are normally distributed with shared 
variance &sigma;<sub>*y*</sub><sup>2</sup>.

## Data Prep

The data for both the response **Y**<sub>*g,t,i*</sub> and the covariate *x*<sub>*ti*</sub> should be preprocessed before 
any analysis. The preprocessing steps are not only crucial, but they are also nontrivial. The preprocessing pipeline provided
is only an outline of a suggestion. Research using other sources, such as [Andrew Jahn's Videos](https://www.youtube.com/channel/UCh9KmApDY_z_Zom3x9xrEQw) 
and [Lindquist et al.'s paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3318970/) on the haemodynamic response 
function are critical. [F. DuBois Bowman's paper](https://www.annualreviews.org/doi/full/10.1146/annurev-statistics-022513-115611) 
is also an excellent resource.

### Response Tensor
When pulling an fMRI scan from a data source, preprocessing is always important in order to make it possible 
to perform any meaningful analyses. As this project is meant to be used with multi-subject data, the suggested 
preprocessing steps are:

  1. Extract the brain image (skullstripping)
  2. Map each subject's brain to some standard (like the Montreal Neurological Institute [MNI] standard)
  3. Correct for motion.
  4. Smooth the image with some kind of statistical method. One viable option is Gaussian smoothing with a full-width half-maximum (FWHM) of 5mm.
  5. Use some brain image software (like FSL) to extract different regions of interest (ROIs) from the brain scans.
  6. Gather the scans for each individual and time, and combine them into a multidimensional array, or *tensor*.
