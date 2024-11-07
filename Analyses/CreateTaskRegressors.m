function CreateTaskRegressors(cfg)
% function CreateTaskRegressors(cfg)
%
% Creates regressors for a GLM in which the behaviour is not modelled, so
% only the condition (congruent vs incongruent) x stimulus (present vs
% absent). Also model responses as nuisance variables. 

nSub = length(cfg.subjects);

for sub = 1:nSub

    outDir = fullfile(cfg.out,cfg.subjects{sub},'Regressors','TaskRegressors');
    if ~exist(outDir,'dir'); mkdir(outDir); end

    dataDir = fullfile(cfg.root,cfg.subjects{sub},'Behaviour');
    dataFiles = str2fullfile(dataDir,'MT*.mat');
    
    if ~iscell(dataFiles); tmp = dataFiles; clear dataFiles; dataFiles{1} = tmp; clear tmp; end
    nRuns     = length(dataFiles);

    for run = 1:nRuns
        saveName = fullfile(outDir,sprintf('run_%d',run));

        % load the info
        load(dataFiles{run},'T','R','trials','blocks','miniblocks','C')

        p_ori = miniblocks; clear miniblocks;
        i_ori = reshape(repmat(blocks,1,length(p_ori)/length(blocks))',length(p_ori),1);
        clear blocks; blocks(:,1) = p_ori; blocks(:,2) = i_ori;
        nBlcks = size(T.presTimes,1); 
        cong = double(blocks(:,1)==blocks(:,2)); % congruency

        %% Nuisance regressors
        Nuisance_names = {'Detection screen','Vividness screen','Incorrect trials'};

        % detection and vividness responses
        t = 1;
        for b = 1:nBlcks
            for m = 1:length(trials)
                % detection
                Nuisance{1}(t,1) = T.presTimes(b,m,2)-T.starttime; % onset screen
                Nuisance{1}(t,2) = R(b,m,4); % response time
                Nuisance{1}(isnan(Nuisance{1}(:,2)),2) = 10; 

                % vividness
                Nuisance{2}(t,1) = T.presTimes(b,m,3)-T.starttime; % onset screen
                Nuisance{2}(t,2) = R(b,m,2); % response time
                Nuisance{2}(isnan(Nuisance{2}(:,2)),2) = 10;

                t = t+1;
            end
        end

        % incorrect imagery blocks 
        ima_check = reshape(repmat(C,1,length(p_ori)/length(C))',length(p_ori),1);
        if sum(ima_check==0)~=0
            inc_blocks = ima_check==0;
            timings    = T.presTimes(inc_blocks,:,1)-T.starttime;
            Nuisance{3}(:,1) = timings(:);
            Nuisance{3}(:,2) = ones(length(timings(:)),1)*2;

            % remove incorrect blocks
            T.presTimes(inc_blocks,:,:) = [];
            cong(inc_blocks) = [];
            trials(inc_blocks,:) = [];
        end

        %% Factors of interest regressors
        factor1 = {'cong','inco'};
        factor2 = {'pres','abs'};
        nCnd    = length(factor1)*length(factor2);

        Onsets  = cell(nCnd,1);
        Stim    = cell(nCnd,1);
        Durs    = cell(nCnd,1);
       
        count = 1;
        for f1 = 1:length(factor1)
            trls = trials(cong==(f1-1),:);
            timings = T.presTimes(cong==(f1-1),:,1)-T.starttime;

            for f2 = 1:length(factor2)
                idx = trls==(f2-1);
                Onsets{count} = timings(idx);
                Durs{count}   = ones(sum(idx(:)),1)*2;
                Stim{count}   = sprintf('%s_%s',factor1{f1},factor2{f2});

                count = count+1;
            end
        end

        % save regressors
        save(saveName,'Onsets','Durs','Stim','Nuisance','Nuisance_names')

        clear trials blocks Onsets Durs Stim Nuisance R T
    end

end