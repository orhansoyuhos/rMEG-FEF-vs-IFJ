%% Project Title: Functional Connectivity Fingerprints of Frontal Eye Field and Inferior Frontal Junction
%%% Skip this part for oPEC and PDC connectivity measures! 
%%% Update Date: 18-March-2022


%% Compatible with fieldtrip-20210411 and MATLAB2020a

%% Output

% outputs per subject connectivity matrices for group level analysis
% AND finds the top-K scouts that have the best connectivity measures with the selected RIOs per subject

%% Inputs
% subjectid = cimec.subjectid; % e.g. '111514'
% top_k = cimec.top_k; % e.g. 5
% trialDuration = cimec.trialDuration; % e.g. 10 (seconds)
%
% functionalConnectivity = cimec.functionalConnectivity; %e.g. '...\HCP_Subjects\Outputs\functional_connectivity\'
%
% % Output file
% outputdir = cimec.outputdir; % '...\HCP_Subjects\Outputs\connectivity_results\';

%% Function

function results_connectivity(cimec)

%% Info about the data
subjectid = cimec.subjectid; % e.g. '111514'
top_k = cimec.top_k; % e.g. 5
trialDuration = cimec.trialDuration; % e.g. 10 (seconds)

functionalConnectivity = cimec.functionalConnectivity; %e.g. '...\HCP_Subjects\Outputs\functional_connectivity\'

% Output file
outputdir = cimec.outputdir; % '...\HCP_Subjects\Outputs\connectivity_results\';

%% Frequency bands
band_names  =   {  'delta' ,   'theta' ,   'alpha' ,   'beta'  ,   'gamma'  };
freq_bands  =   [    1 4   ;     4 8   ;    8 13   ;    13 30  ;    30 100  ];

N_bands = length(band_names);

%% Loop over frequency bands
conn_results = {};

for band = 1:N_bands
    
    %% Load frequency band
    inputfile = fullfile(functionalConnectivity, [subjectid, '_functionalConnectivity_', band_names{band}, sprintf('_%ds.mat', trialDuration)]);
    load(inputfile, 'avgFreq_funConn')
    
    conn_metric = avgFreq_funConn.conn_metric;
    
    %% Organize the pairs of scouts' connectivity value
    conn_results(band).band = band_names{band};
    conn_results(band).name = unique(avgFreq_funConn.labelcmb(:,1));
    
    N = length(conn_results(band).name);
    max_k = size(avgFreq_funConn.labelcmb, 1);
    step = max_k/N;
    k = [0:step:max_k];
    
    conn_results(band).data = zeros(step, N);
    conn_results(band).pair = cell(step, N);
    
    for i = 1:N
        
        if strcmp(conn_metric, 'iCOH')
            % for imaginary part of coherency
            conn_results(band).data(:,i) = [avgFreq_funConn.cohspctrm(k(i)+1:k(i)+step)]';
            
        elseif strcmp(conn_metric, 'dwPLI')
            % https://mailman.science.ru.nl/pipermail/fieldtrip/2011-June/003892.html
            % https://mailman.science.ru.nl/pipermail/fieldtrip/2019-June/039226.html
            % https://mailman.science.ru.nl/pipermail/fieldtrip/2016-February/022995.html
            % https://mailman.science.ru.nl/pipermail/fieldtrip/2014-January/020309.html
            
            % debiased weighted phase lag index
            conn_results(band).data(:,i) = [avgFreq_funConn.wpli_debiasedspctrm(k(i)+1:k(i)+step)]';
            
            
        end
        conn_results(band).pair(:,i) = [avgFreq_funConn.labelcmb(k(i)+1:k(i)+step, 2)]';
        
    end
    
    %% Top K regions that has high connectivity with selected RIOs
    conn_results(band).top_k.name = cell(top_k, N);
    conn_results(band).top_k.data = zeros(top_k, N);
    
    for i = 1:N
        
        values = conn_results(band).data(:, i);
        [max_top, idx_top] = maxk(values, top_k);
        conn_results(band).top_k.data(:,i) = max_top;
        conn_results(band).top_k.name(:,i) = conn_results(band).pair(idx_top, i);
        
    end
    
    conn_results(band).band_names = band_names;
    conn_results(band).freq_bands = freq_bands;
    conn_results(band).label = avgFreq_funConn.label;
end

%% Save the results

% .mat file
mkdir(outputdir);
save(fullfile(outputdir, [subjectid, '_connectivityResults_', sprintf('%ds.mat', trialDuration)]), 'conn_results', '-v7.3');

%% Excel file

% results
for idx_band = 1:length(conn_results(1).band_names)
    
    band_name = conn_results(idx_band).band;
    results(idx_band).band_name = band_name;
    
    for idx_RIO = 1:size(conn_results(idx_band).name, 1)
        
        RIOs_name = conn_results(idx_band).name{idx_RIO};
        top_k = conn_results(idx_band).top_k.name(:, idx_RIO);
        data = conn_results(idx_band).top_k.data(:, idx_RIO);
        
        results(idx_band).top_k{idx_RIO} = top_k;
        results(idx_band).data{idx_RIO} = data;
        results(idx_band).RIOs_name{idx_RIO} = RIOs_name;
    end
end

% make a table
excel_name = fullfile(outputdir, [subjectid, '_topK_', sprintf('%ds.xlsx', trialDuration)]);

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
    writetable(make_table, excel_name, 'Sheet', results(idx_band).band_name, 'Range', 'C3:E40');
    
end

disp('The results are saved.')

end