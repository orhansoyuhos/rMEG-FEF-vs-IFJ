## Analysis_Pipeline folder includes MATLAB code to implement functional connectivity analysis.

AnalysisPipeline_Part1.m script preprocesses the MEG recordings and divides them into desired epoch lengths. It is compatible with megconnectome 3.0, fieldtrip-r10442 and MATLAB2012b.

AnalysisPipeline_Part2.m is an END-to-END script to implement group-level functional connectivity analysis based on the preprocessed MEG data. It is compatible with brainstorm (Version: 28-May-2021), fieldtrip-20210411 and MATLAB2020a.

Before running the code, please fill the directory information in the load_path.m script inside the helper_fun folder to run the code.

### Ground-truth analyses

oPEC: orthogonalized Power Envelope Correlation

![3](https://user-images.githubusercontent.com/44211738/159386267-29a470da-98c1-4d02-bc8b-d00d20b22017.PNG)

