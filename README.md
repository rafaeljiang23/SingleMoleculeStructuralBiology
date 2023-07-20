# SingleMoleculeStructrualBiology
This is a workflow to analyze isolated single membrane transporters (GltPh) from HS-AFM experiments. The codes are developed in BIO-AFM-LAB at Weill Cornell Medicine

Developer: Yining Jiang

Publication: XXX

User should email to the corresponding author of the paper for details about this work: Professor Simon Scheuring (sis2019@med.cornell.edu)

NOTE: at this stage, the code is for peer review only.

## System requirements:
1. Operating system for code development : macOS Big Sur Version 11.7.8
2. Software for code development: MATLAB (MathWorks) 2021b
3. Additional add-ons: MIJI
4. Non-standard hardware: N/A

## Installation instructions: 
1. The codes require installation of MATLAB (MathWorks) 2021b. An installation guide can be found at: https://www.mathworks.com/help/install/.
2. MIJI is recommanded (not required) for visualizing data. An installation guide can be found at: https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab.
3. The installation should take less than one hour on a "normal" desktop computer

## Demo
### Instructions to run on data
These codes comprise of two parts:
#### 1. Sorting protomer observations (PO)
##### Main scripts:
1. PCA_1_create_PO_objects.m
2. PCA_2_locate_PO_COM.m
3. PCA_3_collect_PO_stats.m
4. PCA_4_sort_IFS_OFS.m
5. PCA_5_sort_IFSo_IFSc.m
##### Helper functions:
1. fill_protomers_idx_nan.m
2. energylanscape.m

##### Instruction
User should follow the order to run the main scripts. Operation details are provided within the scripts.

#### 2. Constructing (3D) localization AFM (LAFM) maps
##### Main scripts:
1. tD1_LAFM_find_max.m
2. tD2_LAFM_organize_detections.m
3. tD3_LAFM_core.m

##### Helper functions:
1. tDAFM_voxels.m
2. make_3D_LAFM_kernel1e.m

##### Instruction
User should follow the order to run the main scripts. Operation details are provided within the scripts.

### Test data
A test data (HS-AFM movie of an isolated cytoplasmic GltPh for ~30s, ~300 frames) is provided in 'PO_test.mat file'. Two outcome files are provided ('PO_test_sort_outcome.mat' and 'PO_test_lafm_outcome.mat') to give main outcome files from the test data. Note that the sorting result may slightly vary due to the randomized intial conditions used in the Gaussian Mixture Models clustering operation, but the outcomes should not change much. Besides, the LAFM detection density (normalized) volume data, 'voxels_hsn', in 'PO_test_lafm_outcome.mat' corresponds to the IFS open state of GltPh. 
These tests are expected to run for less than 30 minutes for demo on a "normal" desktop computer following the instructions provided in the main scripts.

## Instruction for use
User should flatten and align the HS-AFM data (of an isolated individual membrane protein), which should be read in MATLAB as the 'data' variable. User should adjust the parameters in the main scripts accourding to the nature of their membrane proteins of interest.
