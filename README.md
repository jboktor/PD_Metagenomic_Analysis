# Metagenomic Analysis of the Parkinson's Disease Microbiome (Under Construction)
### This repository recreates the shotgun metagenomic analysis presented in : paper link

***

##### Badges will go here

## Table of Contents
1. [Background](#Background)
2. [Requirements]("#Requirements")
3. [Workflow](#Workflow)

## Background
Some text here on Parkinsons disease: Add Abstract


## Requirements
All software used for this analysis is open source and freely available to the public. 
The majority of this analysis takes place in R-studio. Certain packages require an R version >= 3.6. 
We recommend updated any R version below this prior to running this analysis to 4.0.1. - "See Things Now"

You may use the links below to download R and R-Studio:

1. Download [R](https://www.r-project.org/) 
2. Download [R-studio](https://rstudio.com/products/rstudio/download/)

In addition, the FlashWeave Probablistic Graphical Models utilze Julia, Python, JupyterNotebook, 
html, CSS, and javascript, follow the intructions below to download the necessary software:

1. Download [python](https://www.python.org/downloads/)
2. Download [JupyterNotebook](https://jupyter.org/install)
3. Download [Julia](https://julialang.org/) (see below)

We recommend downloading the binary from https://julialang.org/downloads/. Julia 1.0 or above are currently supported by FlashWeave.

> To call julia from the command line, add your downloaded path to your .bash_profile:
> Update .bash_profile with the following (__BUT__ Replace quoted section with your own download location/version):

`PATH="/Applications/Julia-1.4.app/Contents/Resources/julia/bin/:${PATH}"
export PATH`

> This line sources your .bashrc file (also add to .bash_profile)

`if [ -f $HOME/.bashrc ]; then
    . $HOME/.bashrc
fi`

## Workflow:
list of R-scripts to run in particular order

1. Set-up/Installation: Run configure.enviornment.R and ensure all packages load sucessfully, if any errors present themselves ensure proper version of R is loaded.
2. create_phyloseq_obj.R
2. PERMANOVA_Analysis.R -> PERMANOVA_Viz.R
3. Beta_Diversity.R
4. Alpha_Diversity.R
5. MaAsLin2_Analysis.R -> MaAsLin2_QC.R
6. Differential_Abundance_Viz_Taxa_Figure_2.R
7. Differential_Abundance_Viz_Functional.R
8. FlashWeave_input_prep.R

## Acknowledgements


## License
A short snippet describing the license (MIT, Apache etc)
MIT License 2020 jboktor

Any questions contact: jboktor@caltech.edu
