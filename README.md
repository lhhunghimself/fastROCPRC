# fastROCPRC
C++ code to calculate area under the curve for ROC and precision recall cures
##Acknowledgements
fastROCPRC based upon R code by Chris Fraley from the [networkBMA package](https://www.bioconductor.org/packages/release/bioc/html/networkBMA.html)
##Summary
Will calculate aurea under the curve (AUC) using the trapezoid rule for Receiver Operator Curves (ROC) and Precision Recall curves (PRC). Also provides tables necessary to plot the curves and summary statistics. The R code is way too slow for evaluating large datasets (hours to days) whereas on the same datasets fastROCPRC completes the calculations in seconds.
