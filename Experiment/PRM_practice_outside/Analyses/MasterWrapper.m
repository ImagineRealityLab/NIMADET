addpath(genpath('Analyses'));
dataDir = 'Results';

dataName = '01';

%% Analyze data
cfg = [];
cfg.nMB = 2; % per block
cfg.dir = dataDir;
cfg.files{1} = dataName;
cfg.plot = true;

analyse_data(cfg);