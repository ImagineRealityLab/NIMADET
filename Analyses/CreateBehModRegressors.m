function CreateBehModRegressors(cfg)
% function CreateBehModRegressors(cfg)
%
% We model the same factors as 'CreateTaskRegressors' but now we
% additionally include the behavioural responses - the reality judgement
% and the vividness rating - as parametric modulators. Because for some
% participants for some runs there is not a lot of variation, we later
% concatenate the runs. Parametric regressors will be orthogonalized with
% respect to the main effect, averaged over runs, not to each other. 
nSubs = length(cfg.subjects);

% column info of the B matrix with behavioural data
cP = 1; cI = 2; cPr = 3; cR = 4; cV = 5; cIC = 6; cTS = 7;
cTR_s = 8; cTR_r = 9; cTV_s = 10; cTV_r = 11;
% Information on columns:
% cP = perception orientation 
% cI = imagery orientation
% cPr = stimulus presence (1=present, 0=absent)
% cR = detection response (reality judgement) (1=present, 0=absent)
% cV = vividness rating
% cIC = imagery check (1=correct, 0 = incorrect)
% cTS = onset stimulus
% cTR_s = onset detection screen
% cTR_r = onset detection response
% cTV_s = onset vividness screen
% cTV_r = onset vividness response 


%% Loop over subjects
for sub = 1:nSubs

    %% Get the behavioural data
    load(fullfile(cfg.root,cfg.subjects{sub},'Behaviour',['RMs_' cfg.subjects{sub} '.mat']),'responseMappings')

    % get the data
    data = str2fullfile(fullfile(cfg.root,cfg.subjects{sub},'Behaviour'),...
        'MT_*.mat');
    nRuns = length(data);

    output_dir = fullfile(cfg.out,cfg.subjects{sub},'Regressors','Behaviour_matrix');
    if ~exist(output_dir,'dir'); mkdir(output_dir); end

    if ~exist(fullfile(output_dir,'run_2.mat'),'file')
        
        for r = 1:nRuns

            load(data{r},'R','C','trials','blocks','miniblocks','T');

            p_ori = miniblocks; clear miniblocks;
            i_ori = reshape(repmat(blocks,1,length(p_ori)/length(blocks))',length(p_ori),1);
            ima_check = reshape(repmat(C,1,length(p_ori)/length(C))',length(p_ori),1);

            nBlcks = size(T.presTimes,1);
            nTrls  = length(trials);

            B      = nan(nBlcks*nTrls,11);

            if responseMappings(r) == 2
                R(:,:,1) = 5-R(:,:,1); % reverse vividness rating
            end

            count = 1;
            for b = 1:nBlcks
                for t = 1:nTrls

                    B(count,cP) = p_ori(b);
                    B(count,cI) = i_ori(b);

                    B(count,cPr) = trials(b,t);

                    B(count,cR) = R(b,t,3);
                    B(count,cV) = R(b,t,1);

                    B(count,cIC) = ima_check(b);

                    B(count,cTS) = T.presTimes(b,t,1)-T.starttime;

                    B(count,cTR_s) = T.presTimes(b,t,2)-T.starttime; % detection screen
                    B(count,cTR_r) = R(b,t,4); % response time

                    B(count,cTV_s) = T.presTimes(b,t,3)-T.starttime; % detection screen
                    B(count,cTV_r) = R(b,t,2); % response time

                    count = count+1;
                end
            end

            % save
            save(fullfile(output_dir,sprintf('run_%d',r)),'B')

            clear R C trials blocks
        end
    end

    %% Make the regressors
    reg_dir = fullfile(fileparts(cfg.root),'Results',cfg.subjects{sub},'Regressors','BehModRegressorsIC');
    if ~exist(reg_dir,'dir'); mkdir(reg_dir);end

    % first calculate global mean per factor for orthogonalization
    behFiles = str2fullfile(output_dir,'run*.mat');
    dat = [];
    for b = 1:length(behFiles)
        load(behFiles{b},'B');
        dat = cat(1,dat,B); clear B
    end
    dat(dat(:,cIC)==0,:) = []; % remove incorrect trials

    % calculate means for mean-centring
    mV = nan(2,2); mR = nan(2,2);
    cong = dat(:,cP)==dat(:,cI);
    for c = 1:2
        for p = 1:2
            idx = cong==(c-1) & dat(:,cPr)==(p-1);
            mV(c,p) = mean(dat(idx,cV));
            mR(c,p) = mean(dat(idx,cR));
        end
    end

    if ~exist(reg_dir,'dir'); mkdir(reg_dir);end

    for r = 1:nRuns

        % get the behavioural info
        load(fullfile(output_dir,sprintf('run_%d',r)),'B')
        IC   = B(:,cIC)==1; 

        % define nuisance regressors
        Nuisance_names{1} = 'Detection rating';
        Nuisance{1}(:,1)  = B(:,cTR_s);
        Nuisance{1}(:,2)  = B(:,cTR_r);

        Nuisance_names{2} = 'Vividness rating';
        Nuisance{2}(:,1)  = B(:,cTV_s);
        Nuisance{2}(:,2)  = B(:,cTV_r);
    
        Nuisance_names{3} = 'Incorrect trials';
        Nuisance{3}(:,1)  = B(~IC,cTS);
        Nuisance{3}(:,2)  = ones(length(B(~IC,cTS)),1).*2;
        
        % remove incorrect trials
        B(~IC,:) = []; 
        
        cong = B(:,cP)==B(:,cI);

        Stim = cell(4,1); Onsets = cell(4,1); Durs = cell(4,1); 
        Viv  = cell(4,1); RJ     = cell(4,1);

        congruency = {'inco','cong'};
        stimulus   = {'abs','pres'};

        % task regressors and pmods
        count=1;
        for c = 1:2
            for p = 1:2
                idx  = cong==(c-1) & B(:,cPr)==(p-1);
                
                Stim{count} = [congruency{c} '_' stimulus{p}];
                Onsets{count} = B(idx,cTS);
                Durs{count}   = ones(sum(idx),1)*2;
                Viv{count}    = B(idx,cV)-mV(c,p); % global mean corrected
                RJ{count}     = B(idx,cR)-mR(c,p); % global mean corrected
                count = count+1;
            end
        end        

        % save 
        save(fullfile(reg_dir,sprintf('run_%d.mat',r)),'Onsets','Durs','Stim','Viv','RJ','Nuisance','Nuisance_names');
        clear Onsets Durs Stim Viv Nuisance Nuisance_names B IC

    end


end
