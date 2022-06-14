# Zebrafish Melanoma

Contains code to analyze fluorescent zebrafish melanoma images and run the GRM model.

Raw image data (amounting to several GBs) as well as certain files required to run the code are available in a public repository (https://doi.org/10.7910/DVN/PQHNU0). This github repository contains the MATLAB code but not the data.

Included in the dataset referenced above are MATLAB code to solve the Growth-Reduction-Metastasis (GRM) model, images of fluorescent melanoma in zebrafish, MATLAB code to analyze the data to obtain tumor size distributions and MATLAB code to find parameters that optimally fit the size distribution data.

The following losely categorizes the content in the dataset based on their purposes:

Raw Image Files:
The files named White Lab\td######_2.tar are compressed images of zebrafish generated from a Zeiss V16 microscope, where the 6 digit number in the name after 'td' is the transplantation date of each batch of fish in ddmmyy format. 

Image Analysis Code:
For example, using Heilmann-fish-image-analysis-mod\CreateFishStructNested_v4_1.m takes the folders containing the raw images from a batch of fish and analyzes using a multistep process where images from the same fish at different time points are aligned, tumors are identified, and the sizes of individual tumors at each time point are collected. The code Heilmann-fish-image-analysis-mod\runCreateFishStruct_v2_4_all.m supplies the right inputs to CreateFishStructNested_v4_1.m to analyze all the image data in the dataset. 
An example of output from this analysis are the files named Heilmann-fish-image-analysis-mod\ct10_lt015_ht04_021720_nosat_B013017NR1_VENTRAL_5e5_summary_*.mat where the B013017NR1 indicates the batch of fish the data was generated from. Aditionally, reviseTumorImageDataHelper_v1_2.m uses Heilmann-fish-image-analysis-mod\reviseTumorImageData_v1_2.m to allow a user to manually review each tumor identified in from the pictures and make some limited revisions such as removing an erroneously identified tumor from the dataset.

Code to generate Tumor Size Distributions from Data
The file makeFishDataDistributions_v1_13 creates distributions of tumor size vs time from the output of the image analysis code. Several different versions of these distributions are constructed using different rules. The directory 'ct10_lt015_ht04_021720_newact_s2_011321' contains example distributions.

Code to Solve the GRM Model
The function named getTransformedDistribution_CTC_v1_1.m is used to solve the model. Subroutines such as in runModelHelper_zf_fits_v2_9.m are used to provide parameters and initial conditions for solving the model and also to create a framework to compare model results to data.

Optimize Parameters
Subroutines such as optimizeWithGetSelectedOptima_v4_4.m optimize parameters to fit the experimental data based on objective functions defined in getFitter_v4_1.m. For a permutation test, the subroutine in optimizePermutation_Stage1_v1_0.m was used to optimize the model against several (permuted) sets of data in parallel.

Code to Create Figures for Publication
Figures and statistical tests for the article, "Mathematical tumor model, zebrafish melanoma validation showing female advantage and new potential dormancy/recurrence mechanism", were generated using paper1_PlotCompilations_v1_4.m
