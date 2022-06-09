%% Project Title: Functional Connectivity Fingerprints of Frontal Eye Field and Inferior Frontal Junction
%%% Functional connectivity results per subject.
%%% Update Date: 18-March-2022

%% Output

% computes the connectivity matrices per subject.


%% Inputs
% subjectid = cimec.subjectid; % e.g. '111514'
% seedRegions = cimec.seedRegions; % e.g. {'L_FEF_ROI L', 'L_IFJa_ROI L', 'L_IFJp_ROI L', 'R_FEF_ROI R', 'R_IFJa_ROI R', 'R_IFJp_ROI R'};
% trialDuration = cimec.trialDuration; % e.g. 10 (seconds)
%
% redefine = cimec.redefine; % e.g. true
%
% conn_metric = cimec.conn_metric; % e.g. iCOH; oPEC; dwPLI
%
% scoutsTimeSeries = cimec.scoutsTimeSeries; % e.g. '...\HCP_Subjects\Outputs\scouts_timeSeries\'
%
% outputdir = cimec.outputdir; % e.g.  '...\HCP_Subjects\Outputs\connectivity_matrices\';

%% Function

function functional_connectivity(cimec)

%% Info about the data
subjectid = cimec.subjectid; % e.g. '111514'
seedRegions = cimec.seedRegions; % e.g. {'L_FEF_ROI L', 'L_IFJa_ROI L', 'L_IFJp_ROI L', 'R_FEF_ROI R', 'R_IFJa_ROI R', 'R_IFJp_ROI R'};
trialDuration = cimec.trialDuration; % e.g. 10 (seconds)

redefine = cimec.redefine; % e.g. true

scoutsTimeSeries = cimec.scoutsTimeSeries; % e.g. '...\HCP_Subjects\Outputs\scouts_timeSeries\'
inputfile = fullfile(scoutsTimeSeries, [subjectid, '_scoutsTimeSeries_', sprintf('%ds.mat', trialDuration)]);

conn_metric = cimec.conn_metric;

% Output file
outputdir = cimec.outputdir; % e.g.  '...\HCP_Subjects\Outputs\functional_connectivity\';

%% Frequencies of interest
% frequency bands
band_names  =   {  'delta' ,   'theta' ,   'alpha' ,   'beta'  ,   'gamma'  };
freq_bands  =   [    1 4   ;     4 8   ;    8 13   ;    13 30  ;    30 100  ];

N_bands = length(band_names);

% step
foi_step    =   [    0.5   ;     0.5   ;     0.5   ;      1    ;      2     ];
% smoothing
tapsmofrq   =   [    0.5   ;     0.5   ;     0.5   ;      1    ;      5     ];

%% Pipeline
%% Load FT_data (scoutsTimeSeries)
load(inputfile, 'FT_data')

% selected scouts
scouts = {};

for r = 1:length(seedRegions)
    for i = 1:length(FT_data.label)
        scouts{end+1, 1} = seedRegions{r};
        scouts{end, 2} = FT_data.label{i};
    end
end

%% Loop over frequency bands
if redefine == true
    cfg         = [];
    cfg.length  = trialDuration;
    cfg.overlap = 0;
    data        = ft_redefinetrial(cfg, FT_data);
else
    data        = FT_data;
end

for i = 1:N_bands
    
    %% power spectral density (PSD)
    cfg             = [];
    %cfg.foilim     = [freq_bands(i,1) freq_bands(i,end)];
    cfg.foi         = [freq_bands(i,1) : foi_step(i) : freq_bands(i,end)];
    cfg.method      = 'mtmfft';
    cfg.output      = 'fourier';
    cfg.tapsmofrq   = tapsmofrq(i);
    freq            = ft_freqanalysis(cfg, data);
    
    %     dir_freq = fullfile(outputdir, '_freqanalysis');
    %     mkdir(dir_freq);
    %     save(fullfile(dir_freq, [subjectid, '_freqanalysis_', band_names{i}, sprintf('_%ds.mat', trialDuration)]), 'freq', '-v7.3');
    
    
    %% connectivity analysis
    % https://www.fieldtriptoolbox.org/reference/ft_connectivityanalysis/

    if strcmp(conn_metric, 'iCOH')
        % imaginary part of coherency
        cfg            = [];
        cfg.method     = 'coh';
        cfg.complex    = 'absimag';
        cfg.channelcmb = scouts;
        conn           = ft_connectivityanalysis(cfg, freq);

    elseif strcmp(conn_metric, 'oPEC')
        % orthogonalized power envelope correlation
        cfg            = [];
        cfg.method     = 'powcorr_ortho';
        conn           = ft_connectivityanalysis(cfg, freq);
        
    elseif strcmp(conn_metric, 'PDC')
        % partial directed coherence
        cfg            = [];
        cfg.method     = 'pdc';
        conn           = ft_connectivityanalysis(cfg, freq);
        
    elseif strcmp(conn_metric, 'dwPLI')
        % https://mailman.science.ru.nl/pipermail/fieldtrip/2011-June/003892.html
        % https://mailman.science.ru.nl/pipermail/fieldtrip/2019-June/039226.html
        % https://mailman.science.ru.nl/pipermail/fieldtrip/2016-February/022995.html
        % https://mailman.science.ru.nl/pipermail/fieldtrip/2014-January/020309.html
        
        % debiased weighted phase lag index
        cfg            = [];
        cfg.method     = 'wpli_debiased';
        cfg.channelcmb = scouts;
        conn           = ft_connectivityanalysis(cfg, freq);
        conn.wpli_debiasedspctrm(conn.wpli_debiasedspctrm<0) = 0;
        
    end

    clear freq;
    
    %% average over frequencies
    cfg = [];
    cfg.frequency = 'all';
    cfg.avgoverfreq = 'yes';
    avgFreq_funConn = ft_selectdata(cfg, conn);
    
    clear conn;
    
    %% save funConn
    avgFreq_funConn.label = FT_data.label;
    avgFreq_funConn.conn_metric = conn_metric;
    save(fullfile(outputdir, [subjectid, '_functionalConnectivity_', band_names{i}, sprintf('_%ds.mat', trialDuration)]), 'avgFreq_funConn', '-v7.3');
    clear avgFreq_funConn
end

end
