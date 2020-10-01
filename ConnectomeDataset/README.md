# Analyses on the Connectome Dataset 

## Getting the data

The publicly available connectome dataset can be downloaded from the [Human Connectome Database](https://db.humanconnectome.org). 

We considered 889 subjects that completed all the tasks. For each subject, we downloaded the data of the MOTOR experiment and the data of the RELATIONAL experiment. The two .mat files `SubjectIDMOTOR.mat` and `SubjectIDRELATIONAL.mat` contain the IDs of the subjects we considered in the analysis (note: the subjects are the same in both analyses).

## Importing the data

The scripts `preparesurfacedataMOTOR.m` and `preparesurfacedataRELATIONAL` assume that you have downloaded the data of these subjects and placed the compressed files in a folder (assigned to the variable `dirdata` in the script). To extract the relevant .nii files, we used [7-zip](https://www.7-zip.org/) (freely available). In order to import the .nii files, a [cifti reader](https://github.com/Washington-University/cifti-matlab) has to be added to the current path of Matlab. The imported data are saved in a large .mat file (`dataConnectomeMotor.mat` and `dataConnectomeRelational.mat`) not present in this repository (~ 1GB). These files can be retrieved by downloading the 889 subjects and running the two scripts.

## Decoding analysis

The script that runs the analyses as described in the paper is `comparexvalconnectome.m` and `comparexvalconnectome_iteration.m`, that requires the surface data of the 889 subjects in .mat format (`dataConnectomeMotor.mat` and `dataConnectomeRelational.mat`). The results are stored in `ResultsMOTOR_1000perm_100iter.mat` and `ResultsRELATIONAL_1000perm_100iter.mat`.  *These files are present in the repository*.

## Displaying the results

The figures displayed in the manuscript can be generated using the scripts `EvaluateResultsConnectomeMOTOR.m` and `EvaluateResultsConnectomeRELATIONAL.m`

