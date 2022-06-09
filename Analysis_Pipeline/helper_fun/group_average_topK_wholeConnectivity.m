%% Project Title: Functional Connectivity Fingerprints of Frontal Eye Field and Inferior Frontal Junction
%%% Group average for oPEC and PDC metrics.
%%% Update Date: 18-March-2022


%% Compatible with fieldtrip-20210411 and MATLAB2020a

function group_average_topK_wholeConnectivity(cimec)

%% Info about the data
subjectids = cimec.subjectids; % e.g. subjectids = {'111514', '140117', '162026', '164636', '191033'};
top_k = cimec.top_k; % e.g. 5
trialDuration = cimec.trialDuration; % e.g. 10 (seconds)

functionalConnectivity = cimec.functionalConnectivity; %e.g. '...\HCP_Subjects\Outputs\functional_connectivity\';

% Output file
outputdir = cimec.outputdir; % '...\HCP_Subjects\Outputs\group_results\';

%% Frequency bands
band_names  =   {  'delta' ,   'theta' ,   'alpha' ,   'beta'  ,   'gamma'  };
freq_bands  =   [    1 4   ;     4 8   ;    8 13   ;    13 30  ;    30 100  ];

N_bands = length(band_names);

%% Loop over frequency bands and subjects
group_results = {};

for band = 1:N_bands
    
    group_results(band).band = band_names{band};
    
    for s = 1:length(subjectids)
        load(fullfile(functionalConnectivity, subjectids{s}, [subjectids{s}, '_functionalConnectivity_', band_names{band}, sprintf('_%ds.mat', trialDuration)]));
        conn_metric = avgFreq_funConn.conn_metric;
        
        if strcmp(conn_metric, 'oPEC')
            % orthogonalized power envelope correlation
            tmp_data = avgFreq_funConn.powcorrspctrm;
            
        elseif strcmp(conn_metric, 'PDC')
            % partial directed coherence
            tmp_data = avgFreq_funConn.pdcspctrm;
        end
        
        group_results(band).data(:,:,s) = tmp_data;
        
        clear tmp_data;
    end
    
    group_results(band).band_names = band_names;
    group_results(band).freq_bands = freq_bands;
    group_results(band).label = avgFreq_funConn.label;
    % mean
    group_results(band).data_perSubject = group_results(band).data;
    group_results(band).data = mean(group_results(band).data, 3);
    
    
    group_results(band).subjectids = subjectids;
end

%% Save the results

% .mat file
id = char(string(length(subjectids)));
mkdir(outputdir);
save(fullfile(outputdir, [id, 'Subjects', '_groupResults_', sprintf('%ds.mat', trialDuration)]), 'group_results', '-v7.3');


disp('The results are saved.')

end