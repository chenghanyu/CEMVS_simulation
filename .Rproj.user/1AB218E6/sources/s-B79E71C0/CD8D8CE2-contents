Complex-Valued EM Variable Selection Algorithm 
================================================

The R scripts implement the complex-valued EM variable selection algorithm (C-EMVS) to reproduce Figure 9 and Figure 10 shown in Section 3.2 of the manuscript. Two different parts of programs are included. The folder **circularEMVS** contains R functions for the circular EMVS algorithm. The associated **circularEMVS_sim2.R** script runs the algorithm and reproduce the figures in the manuscript. The R script **EMVSmainfcn.R** includes the main circular EM algorithm and other functions that are used for analysis and picturing. The general non-circular EMVS functions are in the another folder **CVEMVS** which use the Rcpp package to embed C++ code into the R functions. The associated **CVEMVS_sim2.R** script runs the general non-circular algorithm and generates the activation and strengh maps from CV-fMRI data. The two versions of code lead to the same result, showing that the circular EMVS is the special case of the non-circular EMVS when the relation parameter is zero. **EMVSmanfcn** folder includes all other functions used to generate other figures and for analysis in the paper. 


Notations
------------------------------------------------
In the manuscript, $\gamma$ is used to denote the regression coefficients and $\psi$ for the indicator variables. However, the code follow a convention that uses $\beta$ to denote the coefficients and $\gamma$ is for the indicator variables.


Data
------------------------------------------------
All data sets and functions for data manipulation are in the folder **data**. The original simulated data set *GrandPrelimSimData.mat* is stored in .mat format. The data have tansformed into a text file and RData file that can be read and operated in R.


Figures
------------------------------------------------
The figures generated from **circularEMVS_sim2.R** are saved in the **figures** folder.








