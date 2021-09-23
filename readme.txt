This project contains analysis codes and pre-processed datasets required to replicate the results shown in the article
"A compositional letter code in high-level visual cortex explains how we read jumbled words"

Briefly
1)	The files starting with figure_* contains code necessary to generate all the plots in that figure panel.
2)	Data files used by these codes are stored in L2*.mat files
3)	Some codes require additional library (provided inside lib folder). Add this folder in the MATLAB path.
5)	Folder "exampleROIdef" contains all files necessary to create ROI definitions for an example subject

Detailed
Data
1)	dfmri_LDC.mat 						= All pair-wise dissimilarity (74C2) for each subject across all the 5 ROIs
2)	dis_letter.mat 						= Single letter dissimilarity for upper case and lower case English letters
3)	figure4_ldt_workspace 				= MATLAB Workspace of figure4_ldt
4)	figure5_semanticfeat 				= semantic feature vector for each of the 32 words obtained using word2vec
5)	figureS5_S6_compoundwords_workspace = MATLAB Workspace of figureS5_S6_compoundwords 
6)	figureS7_S9_neural2ISI_workspace 	= MATLAB Workspace of figureS7_S9_neural2ISI 
7)	figureS14_scr_workspace 			= MATLAB Workspace of figureS14_scr 

8)	L2_4letter 		= behavioural data structure for the 4 letter optimized searches (Experiment S4)
9)	L2_6letter 		= behavioural data structure for the 6 letter optimized searches (Experiment S5)
10)	L2_bigram 		= behavioural data structure for the bigram searches (Experiment 2)
11)	L2_ldt 			= behavioural data structure for the Lexical Decision Task (Experiment 5)
12)	L2_letter 		= behavioural data structure for the single letter (A-Z, a-z, and 0-9) searches (Experiment 1)
13)	L2_Rbigram 		= behavioural data structure for the bigram searches (Experiment S2)
14)	L2_scrwrd 		= behavioural data structure for the jumbled word task (Experiment S7)
15)	L2_string 		= behavioural data structure for the multiple string length searches (Experiment S3)
16)	L2_trigram 		= behavioural data structure for the Trigram searches (Experiment S1)
17)	L2_upinv_bigram = behavioural data structure for the upright-inverted Bigram searches (Experiment 3)
18)	L2cb_3letter 	= behavioural data structure for the Trigram used in compound word searches (Experiment 4)
19)	L2cb_6letter 	= behavioural data structure for the compound word searches (Experiment 4)
20)	L2fmri_LDTw 	=  fMRI data structure for the lexical decision task (Experiment 6)

21)	median_droi =  Median (across subjects) pair-wise dissimilarity values across all stimuli estimated in search light analysis
22) NeuronTuning = Position weights for each of the 10 neurons estimated using the compound word experiments 
23) TableS1_unequalstring_neural_workspace = MATLAB Workspace of Table S1 
24) W_gendiss_neural = data required to predict the dissimilarity between any two strings using neural model

Aditional Codes	
1) calculate_medianD = code required to calculate the Median (across subjects) pair-wise dissimilarity values across all stimuli
2) gendiss = Estimate pair-wise dissimilarity between any number of strings using part-sum model
3) gendiss_neural = Estimate pair-wise dissimilarity between any number of strings using neural letter model
4) getvoxind = Code to obtain voxel ids corresponding to each ROI and subject
5) batonfun = function to estimate dissimilarities using SID model
6) varorder = function to rearrange the model parameters obtained from part-sum model to be compatible with ISI/SID model
7) neuralmodel = function to estimate dissimilarities using letter model
8) partmodel = function to estimate dissimilarities using part-model, ISI, SID model (figure S11)
9) partmodel_isi = function to estimate dissimilarities using letter-model, and ISI model (figure S7 and S9)
10) partmodel_neu= function to estimate dissimilarities using letter-model (Table S1)
11) searchlight_rsa_gendis = code to generate search-light data for dissimilarity analysis 
12) searchlight_rsa_comb = code to generate search light plots after loading the dissimilarity values from median_droi
13) Stimulus_neuralActivity = predicting the neural response for a given input string
14) LDC_distance_measurement = This code is used to generate LDC dissimilarities using RSA toolbox and preprocessed data. 