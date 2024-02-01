---
layout: default
permalink: /software/dosresmeta
show_in_nav: false
navigation_weight: 7
---

The dosresmeta R package: Performing multivariate dose-response meta-analysis
========

The `dosresmeta` package is a free and open-source add-on for conducting (multivariate) dose-response meta-analysis with the statistical software environment [R](https://www.r-project.org/).
You can install the package (or upgrade to the newest version) from [CRAN](https://cran.r-project.org/web/packages/dosresmeta/index.html) directly within R by typing:

`install.packages("dosresmeta")`

Extensions to meta-regression models and alternative estimation procedures are included in an update version still under developement, which is avaiable on [GitHub](https://github.com/alecri/dosresmeta) by typing:

`install.packages("devtools")`

`devtools::install_github("alecri/dosresmeta")`

&nbsp;

## Reproducible code

* [*GitHub repository* for **Crippa, A.**, Discacciati A., Bottai M., Spiegelman D., & Orsini N.  “One-stage dose-response meta-analysis for aggregated data 2018. Statistical Methods in Medical Research](https://github.com/alecri/one-stage-dosresmeta)

* [*GitHub repository* for **Crippa, A.**, & Orsini N.  “Dose-response meta-analysis of differences in means" 2016. BMC Medical Research Methodology](https://github.com/alecri/differences-in-mean)

* [*GitHub repository* for Discacciati, A., **Crippa, A.**, & Orsini N.  “Goodness-of-fit tools for dose-response meta-analysis of binary outcomes" 2015. JRSM](https://github.com/anddis/goodness-of-fit-meta-analysis)


&nbsp;

## Analysis Examples

The links below demonstrate how the models, methods, and techniques described in the respective articles/chapters can be applied via the `dosresmeta` package. The examples are run using dosresmeta version 2.0.0 avaiable on GitHub.

* [Orsini, N. Weighted mixed-effects dose-response models for tables of correlated contrasts. Stata Journal. 2021, Vol.21 (2), p.320-347.](https://rpubs.com/alecri/drmeta)

* [Berlin JA, Longnecker MP, Greenland S. Meta-analysis of epidemiologic dose-response data. Epidemiology. 1993 May 1:218-28](http://rpubs.com/alecri/berlin)

* [Orsini, N., Bellocco, R., & Greenland, S. Generalized least squares for trend estimation of summarized dose-response data. 2006 Stata Journal, 6(1), 40.](http://rpubs.com/alecri/glst)

* [Liu, Qin, et al. "A two-stage hierarchical regression model for meta-analysis of epidemiologic nonlinear dose–response data." Computational Statistics & Data Analysis 53.12 (2009): 4157-4167.](http://rpubs.com/alecri/qinliu)

&nbsp;


**Additional (useful) code**

* [Missing number of cases/n. Assuming independence](/downloads/codes/missing cases n.txt)  
* [Loop for selecting best fractional polynomial](/downloads/codes/fp_splines_example.txt)  
* [Leave-one-out dose-response meta-analysis](/downloads/codes/leave1out.txt)  
* [Cumulative dose-response meta-analysis](/downloads/codes/cumulative.txt)  
* [Sensitivity analysis to location of knots in a spline model](/downloads/codes/knots_location.txt)


**Older R script in txt format based on dosresmeta version 1.3.1 available on CRAN.**

* [Coffee consumption and all causes mortality](/downloads/codes/coffee_mort.txt)
* [Alcohol intake and colorectal cancer](/downloads/codes/alcohol_crc.txt)
* [Alcohol intake and lung cancer](/downloads/codes/alcohol_lc.txt)
* [Body mass index and renal cancer](/downloads/codes/ts_glst_bmi_ex.txt)
* [Alcohol intake and cardiovascular risk](/downloads/codes/ts_glst_alcohol_ex.txt)
