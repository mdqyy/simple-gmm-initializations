


This folder contains the evaluations (single and summarized rankings) as described in the paper.


1. Overview of evaluations in this folder
==============================================================================================


summarized-rankings
----------------------------------------------------------------------------------

- Files
	Mappe_k10.xlsx  and
	Mappe-k10-noisy.xlsx
  can be opened with Excel/LibreOffice.
  Contain evaluations with respect to data without noisy and noisy data, respectively.

- Each file contains multiple sheets:

  * overview	
	contains the summarized rankings that we refer to in the paper.
	Each of the rankings is with respect to the average negative-log-likelihood.
	There are 7 summarized rankings, where we added the rankings wrt data sets that are/have
	(1) in the Spherical test sets, 
	(2) in the Elliptical test sets,
	(3) in the Elliptical-Difficult test sets,
	(4) in some test set (i.e., we summed all existing rankings here),
	(5) a separation of 0.5,
	(6) a separation of 1,
	(7) a separation of 2.	

  * spherical-mean
	contains the summarized rankings for each of the Spherical test sets
	with respect to the average negative log-likelihood.
	At the bottom of the sheet there are different summaries (of all data sets and data sets with specified separation).
  * spherical-variance
  	contains the summarized rankings for each of the Spherical test sets
	with respect to the variance of the negative log-likelihood.
	At the bottom of the sheet there are different summaries (of all data sets and data sets with specified separation).
  * elliptical-mean
  	contains the summarized rankings for each of the Elliptical test sets
	with respect to the average negative log-likelihood.
	At the bottom of the sheet there are different summaries (of all data sets and data sets with specified separation).
  * elliptical-variance
  	contains the summarized rankings for each of the Elliptical test sets
	with respect to the variance of the negative log-likelihood.
	At the bottom of the sheet there are different summaries (of all data sets and data sets with specified separation).
  * elliptical-difficult-mean
  	contains the summarized rankings for each of the Elliptical-Difficult test sets
	with respect to the average negative log-likelihood.
	At the bottom of the sheet there are different summaries (of all data sets and data sets with specified separation).
  * ellitpical-difficult-variance
  	contains the summarized rankings for each of the Elliptical-Difficult test sets
	with respect to the variance of the negative log-likelihood.
	At the bottom of the sheet there are different summaries (of all data sets and data sets with specified separation).

- For an explanation on columns C-N, which contain the parameters of a data set, see "Generated Data.txt".
- For an explanation on the abbrevations used for the initialization methods, see "2. Generation Information".

  

single-rankings
----------------------------------------------------------------------------------

- Files
	summary_statistics-k10-spherical-2013-11-19.csv,
	summary_statistics-k10-elliptical-2013-11-19.csv,
	summary_statistics-k10-elliptical-difficult-2013-11-19.csv,
	summary_statistics-noisy-k10-spherical-2013-11-27.csv,
	summary_statistics-noisy-k10-elliptical-2013-11-27.csv, and
	summary_statistics-noisy-k10-elliptical-difficult-2013-11-27.csv
  have csv format.
  Contain the summarized rankings per test set.

- For an explanation on columns C-N, which contain the parameters of a data set, see "Generated Data.txt".
- The remaining columns contain the following information
	for the initial solution ("init_") and the final solution ("final_"),
	i.e., the minimum ("min") and maximum ("max") value, the lower ("q25") and upper ("q75") quartile,
	the variance ("var"), the mean ("mean") and the median ("med").
- For an explanation on the abbrevations used for the initialization methods, see "2. Generation Information".


2. General Information
==============================================================================================

Abbrevations used for the Initialization methods
----------------------------------------------------------------------------------
- Uniform	= Unif.2GMM, UniformMeans2GMM
- K-Means++	= KM++2GMM, AdaptiveMeans2GMM
- Adaptive(1)	= Adapt, Adapt.(1), Alt_AdaptMean+Means2GMM_i, Alt_AdaptMean+Means2GMM_f1
- Adaptive(0.5)	= Adapt.(0.5), Alt_AdaptMean+Means2GMM_f0.5
- Adaptive(0.7)	= Adapt.(0.7), Alt_AdaptMean+Means2GMM_f0.7
- Adaptive(0.9)	= Adapt.(0.9), Alt_AdaptMean+Means2GMM_f0.9
- Agglo.(0.1)	= AgglomerativeInit__sample0.1
- Gonzlaez	= Gonz.2GMM, Gonzalez2GMM
- GonzalezForGMM(0.1)	= Gonz.ForGMM(0.1), GonzalezForGMM_.*_2gmm1_sample0.1
- KwedlosGonzalez(0.1)	= Kwedlo(0.1), GonzalezKwedlo_.*_sample0.1_