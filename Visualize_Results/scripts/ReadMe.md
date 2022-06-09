## The code to visualize the results for predominant functional connectivities.

- The helper_fun folder includes custom functions used in the scripts below. Please fill the directory information inside load_path.m file.

### Main scripts

- RUNME_visualizeResults_functionConn.m: To apply statistical analysis and visualize FEF and IFJ's predominant functional connectivities in the ipsilateral or contralateral sides. It produces the results on fsaverage cortex and circular graphs.  

- RUNME_visualizeResults_functionConn_collapseTargets.m: Same as the above script but irrespective of laterilazation. 

- RUNME_visualizeResults_effectiveConn.m: To apply statistical analysis and visualize FEF and IFJ's direction of interaction with the visual streams or the rest of the brain.

### Mass production of the figures

- helper_saveFig_effectiveConn.m and helper_saveFig_functionalConn.m is for the mass production of the figures outputted from the main scripts above. It helps to save figures across frequency bands and lateralization.
