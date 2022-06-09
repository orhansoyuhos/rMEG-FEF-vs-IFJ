%% Project Title: Functional Connectivity Fingerprints of Frontal Eye Field and Inferior Frontal Junction
%%% Group average for iCOH and dwPLI metrics.
%%% Update Date: 18-March-2022


%% Compatible with fieldtrip-20210411 and MATLAB2020a

function group_average_topK_forSelectedROIs(cimec)

%% Info about the data
subjectids = cimec.subjectids; % e.g. subjectids = {'111514', ...};
top_k = cimec.top_k; % e.g. 5
trialDuration = cimec.trialDuration; % e.g. 10 (seconds)

connectivityResults = cimec.connectivityResults; %e.g. '...\HCP_Subjects\Outputs\connectivity_results\';

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
        load(fullfile(connectivityResults, [subjectids{s}, '_connectivityResults_', sprintf('%ds.mat', trialDuration)]), ...
            'conn_results');
        
        tmp_data = conn_results(band).data;
        group_results(band).data(:,:,s) = tmp_data;
        
        clear tmp_data;
    end
    
    group_results(band).name = conn_results(band).name;
    group_results(band).pair = conn_results(band).pair;
    group_results(band).band_names = conn_results(band).band_names;
    group_results(band).freq_bands = conn_results(band).freq_bands;
    group_results(band).label = conn_results(band).label;
    % mean
    group_results(band).data_perSubject = group_results(band).data;
    group_results(band).data = mean(group_results(band).data, 3);
    
    
    %% Top K regions that has high connectivity with selected RIOs
    N = length(group_results(band).name);

    group_results(band).top_k.name = cell(top_k, N);
    group_results(band).top_k.data = zeros(top_k, N);
    
    for i = 1:N
        
        values = group_results(band).data(:, i);
        [max_top, idx_top] = maxk(values, top_k);
        group_results(band).top_k.data(:,i) = max_top;
        group_results(band).top_k.name(:,i) = group_results(band).pair(idx_top, i);
        
    end 
    
    group_results(band).subjectids = subjectids;
end

%% Save the results

% .mat file
id = char(string(length(subjectids)));
mkdir(outputdir);
save(fullfile(outputdir, [id, 'Subjects', '_groupResults_', sprintf('%ds.mat', trialDuration)]), 'group_results', '-v7.3');


%% Excel file

% results
for idx_band = 1:length(group_results(1).band_names)
    
    band_name = group_results(idx_band).band;
    results(idx_band).band_name = band_name;
    
    for idx_RIO = 1:size(group_results(idx_band).name, 1)
        
        RIOs_name = group_results(idx_band).name{idx_RIO};
        top_k = group_results(idx_band).top_k.name(:, idx_RIO);
        data = group_results(idx_band).top_k.data(:, idx_RIO);
        
        results(idx_band).top_k{idx_RIO} = top_k;
        results(idx_band).data{idx_RIO} = data;
        results(idx_band).RIOs_name{idx_RIO} = RIOs_name;
    end
end

% make a table
excel_name = fullfile(outputdir, [id, 'Subjects', '_topK_', sprintf('%ds.xlsx', trialDuration)]);

for idx_band = 1:length(results)
    
    tmp_band = results(idx_band);
    
    for idx_RIO = 1:length(tmp_band.RIOs_name)
        
        Pair = [tmp_band.top_k{idx_RIO}; {'-------------'}];
        Icoh = [tmp_band.data{idx_RIO}; NaN];
        RIO  = [repmat(tmp_band.RIOs_name(idx_RIO), length(Pair)-1, 1); {'-------------'}];
        
        if idx_RIO == 1
            make_table = table(RIO, Pair, Icoh);
        else
            T = table(RIO, Pair, Icoh);
            make_table = [make_table; T];
        end
    end
    
    % save as an excel file
    writetable(make_table, excel_name, 'Sheet', results(idx_band).band_name, 'Range', 'C3:E200');
    
end

disp('The results are saved.')

end