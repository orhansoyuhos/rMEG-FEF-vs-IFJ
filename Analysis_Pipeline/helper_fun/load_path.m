%% Project Title: Functional Connectivity Fingerprints of Frontal Eye Field and Inferior Frontal Junction
%%% Loading the paths. 
%%% Update Date: 18-March-2022

function path = load_path(directory)

% libraries and code
if strcmp(directory, 'brainstorm')
    path = '/home/orhan.soyuhos/Desktop/Orhan/MATLAB_lib/brainstorm3';
elseif strcmp(directory, 'fieldtrip-r10442')
    path = '/home/orhan.soyuhos/Desktop/Orhan/MATLAB_lib/fieldtrip-r10442'; 
elseif strcmp(directory, 'fieldtrip-20210411')
    path = '/home/orhan.soyuhos/Desktop/Orhan/MATLAB_lib/fieldtrip-20210411'; 
    
elseif strcmp(directory, 'megconnectome-30')
    path = '/home/orhan.soyuhos/Desktop/Orhan/MATLAB_lib/megconnectome-3.0'; 
elseif strcmp(directory, 'workingDir_code')
    path = '/home/orhan.soyuhos/Desktop/Thesis_OrhanSoyuhos/Analysis_Pipeline';
 
% raw data    
elseif strcmp(directory, 'dir_subjects')
    path = '/home/orhan.soyuhos/Desktop/Orhan/HCP_Subjects';
    
% middle and final outputs
elseif strcmp(directory, 'scoutsTimeSeries')
    path = '/home/orhan.soyuhos/Desktop/Orhan/HCP_Subjects/Outputs/scouts_timeSeries/';
elseif strcmp(directory, 'functionalConnectivity')
    path = '/home/orhan.soyuhos/Desktop/Orhan/HCP_Subjects/Outputs/functional_connectivity/';
elseif strcmp(directory, 'connectivityResults')
    path = '/home/orhan.soyuhos/Desktop/Orhan/HCP_Subjects/Outputs/connectivity_results/';
elseif strcmp(directory, 'groupResults')
    path = '/home/orhan.soyuhos/Desktop/Orhan/HCP_Subjects/Outputs/group_results/';

   
end
