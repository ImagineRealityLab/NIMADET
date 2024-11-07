function VividnessResponseDecodingROIs2(cfg)
% trains and tests LDA classifiers with n-fold cross-validation on the task
% data - train on response and test vividness and vice versa

rng(1); % for reproduceability

cfg.nPerm = 25; % per participants
cfg.nBtstrp = 10000; % to create group null distribution

[~,atlas] = read_nii(cfg.atlas);
mask      = atlas>0; atlas = atlas(mask);
nSub      = length(cfg.subjects);

trial_type{1} = 'congruency == 1';
trial_type{2} = 'congruency == 0';
nTT = length(trial_type);

acc   = nan(nSub,nTT);
pAcc = nan(nSub,nTT,cfg.nPerm);

trl_distr = nan(2,nSub,2,2); % before/after, cong, viv, RJ

posterior = nan(nSub,nTT,2);

for sub = 1:nSub

    fprintf('Processing subject %s \n',  cfg.subjects{sub})

    % out dir
    outDir = fullfile(cfg.root,'Results',cfg.subjects{sub},cfg.outDir);
    if ~exist(outDir,'dir'); mkdir(outDir); end

    %check if already done
    if ~exist(fullfile(outDir,'ROI_accuracy.mat'),'file')

    %% prepare the data

    % task data
    dataDir = fullfile(cfg.root,'Results',cfg.subjects{sub},cfg.dataDir);
    nRuns   = length(dir(dataDir))-2;
    D = []; I = []; P = []; R = []; V = []; data = []; run_idx = [];
    for run = 1:nRuns

        fprintf('Reading in task betas for run %d... \n',run)

        beta_niftis = str2fullfile(fullfile(dataDir,sprintf('RUN_%02d',run)),'trl_*.nii');
        nifti_names = cell(length(beta_niftis),1); % re-sort
        beta_names  = cell(length(nifti_names),1);
        for n = 1:length(beta_niftis)
            nifti_names{n} = str2fullfile(fullfile(dataDir,sprintf('RUN_%02d',run)),sprintf('trl_%d_*.nii',n));
            [~,beta_names{n}] = fileparts(nifti_names{n});
        end

        nTrls = length(beta_names);

        % remove blocks failed imagery check
        if cfg.ima_check
            load(fullfile(cfg.root,'Data',cfg.subjects{sub},...
                'Behaviour',sprintf('MT_%s_run%d.mat',cfg.subjects{sub},run)),'C')
            if any(C==0)
                task_ridx = [];
                for c = 1:length(C)
                    idx = (c-1)*(nTrls/length(C))+1:c*(nTrls/length(C));
                    if C(c)==0
                        task_ridx = [task_ridx,idx];
                    end
                end
                beta_names(task_ridx) = [];
                nifti_names(task_ridx) = [];
            end
        end

        nTrls = length(beta_names);
        betas = nan(nTrls,sum(mask(:)));
        for t = 1:nTrls
            [~,Y] = read_nii(nifti_names{t});
            betas(t,:) = Y(mask); clear Y;

            num_idx = [strfind(beta_names{t},'Det_') strfind(beta_names{t},'_Ima')];
            D = [D str2double(beta_names{t}(num_idx(1)+4:num_idx(2)-1))];

            num_idx = [strfind(beta_names{t},'Ima_') strfind(beta_names{t},'_P')];
            I = [I str2double(beta_names{t}(num_idx(1)+4:num_idx(2)-1))];

            num_idx = strfind(beta_names{t},'P_');
            P = [P str2double(beta_names{t}(num_idx+2))];

            num_idx = strfind(beta_names{t},'R_');
            R = [R str2double(beta_names{t}(num_idx+2))];

            num_idx = strfind(beta_names{t},'V');
            V = [V str2double(beta_names{t}(num_idx+1))];
        end
        run_idx = [run_idx; ones(nTrls,1)*run];
        clear SPM beta_names nifti_names
        data = cat(1,data,betas); clear betas
    end

    % define congruency
    congruency = D == I;
    task_oris = unique(D);

    % loop over trial types
    task_data = cell(nTT,1); R_task = cell(nTT,1); % responses
    task_ridx = cell(nTT,1);
    for tt = 1:nTT

        trial_type_idx = eval(trial_type{tt});

        task_data{tt} = data(trial_type_idx,:);
        task_ridx{tt} = run_idx(trial_type_idx);
        R_task{tt}    = [V(trial_type_idx)' R(trial_type_idx)'];        

        % first z-score per run
        X = []; Res = []; RN = [];
        RI = unique(task_ridx{tt}); nRuns = length(RI);
        for r= 1:nRuns

            run = RI(r);

            xrun = task_data{tt}(task_ridx{tt}==run,:);
            rrun = R_task{tt}(task_ridx{tt}==run,:);

            RN = cat(1,RN,ones(length(rrun),1)*run);
            Res  = cat(1,Res,rrun);

            % mean centre per run
            X = cat(1,X,xrun-nanmean(xrun,2));
        end
        task_ridx{tt} = RN;
        task_data{tt} = X;
        R_task{tt} = Res;
    end

    %% Do the decoding per ROI per trial type
    acc_sub       = nan(nTT,1); 
    pAcc_sub      = nan(nTT,cfg.nPerm);
    posterior_sub = nan(nTT,2);

    % loop over trial types
    for tt = 1:nTT

        x = task_data{tt};
        y = R_task{tt}; m = median(y(:,1)); 
        % seperate into low vs high
        viv_cat = y(:,1)>=m;%
        run_idx = task_ridx{tt};

        % save the distribution before balancing
        for v = 1:2
            trl_distr(1,sub,tt,v) = sum(viv_cat==(v-1) & y(:,2)==1)./sum(viv_cat==(v-1));
        end

        % remove nans
        nan_idx = any(isnan(x));
        x(:,nan_idx) = [];

        % remove mean over voxels <<<
        x = x-mean(x,2);

        % classify
        if ~isempty(x) && (length(unique(y(:,1)))>1) && (length(unique(y(:,2)))>1)            

            % cross-decode with completely balanced training set
            Y1 = []; Y2 = []; Ridx = [];
            for r = 1:nRuns
                idx = run_idx == r; ridx = run_idx(idx);
                y1  = viv_cat(idx); y2 = y(idx,2);
                y3  = nan(length(y1),1);
                y3(y1==0&y2==0) = 1;
                y3(y1==0&y2==1) = 2;
                y3(y1==1&y2==0) = 3;
                y3(y1==1&y2==1) = 4;

                if length(unique(y3))==4 % only if we have all classes
                    idx = balance_trials(y3,'downsample');
                    Y1  = cat(1,Y1,y1(cell2mat(idx)));
                    Y2  = cat(1,Y2,y2(cell2mat(idx)));
                    Ridx = cat(1,Ridx,ridx(cell2mat(idx)));
                end
            end

            % save the distribution after balancing
            for v = 1:2
                trl_distr(2,sub,tt,v) = sum(Y1==(v-1) & Y2==1)./sum(Y1==(v-1));
            end

            % cross-validate based on runs
            if ~isempty(Y1)
                yhat     = nan(length(Y2),1);
                yhatP    = nan(length(Y2),cfg.nPerm);

                for rn = 1:nRuns

                    test_idx  = find(Ridx==rn);
                    train_idx = setdiff(1:length(Y1),test_idx);

                    trainx = x(train_idx,:);
                    trainy = Y1(train_idx);

                    testx  = x(test_idx,:);

                    % train vividness - test response
                    decoder = train_LDA(cfg,trainy,trainx');
                    yhat(test_idx) = decode_LDA(cfg,decoder,testx');

                    for p = 1:cfg.nPerm

                        % train vividness - test response
                        permY = trainy(randperm(length(trainy)));
                        decoder = train_LDA(cfg,permY,trainx');
                        yhatP(test_idx,p) = decode_LDA(cfg,decoder,testx');
                    end
                end

                acc_sub(tt) = mean((yhat>0)==Y2);
                pAcc_sub(tt,:) = mean((yhatP>0)==repmat(Y2,1,cfg.nPerm));

                posterior_sub(tt,1) = mean(yhat(Y2==0));
                posterior_sub(tt,2) = mean(yhat(Y2==1));

            end        
           

        end

    end

   % write the results
   if ~exist(outDir,'dir'); mkdir(outDir); end
   save(fullfile(outDir,'ROI_accuracy'),'acc_sub','pAcc_sub','cfg','posterior_sub')

   else
      fprintf('\t already done, loading results... \n')
      load(fullfile(outDir,'ROI_accuracy'),'acc_sub','pAcc_sub','posterior_sub')
   end

    acc(sub,:) = acc_sub;
    pAcc(sub,:,:) = pAcc_sub;
    posterior(sub,:,:) = posterior_sub;
end

% save group-results
group_dir = fullfile(cfg.root,'Results','GroupResults',cfg.outDir);
if ~exist(group_dir,'dir'); mkdir(group_dir); end
save(fullfile(group_dir,'ROI_decoding'),'acc','pAcc','posterior');

% remove nans
nan_idx = squeeze(any(any(isnan(posterior),2),3));%any(isnan(acc'))';

acc(nan_idx,:) = []; pAcc(nan_idx,:,:) = [];
posterior(nan_idx,:,:) = [];

%% Plot the results

figure;
subplot(2,1,1)
barwitherr(squeeze(std(acc))./sqrt(nSub),squeeze(mean(acc)));
ylim([0.45 0.6]); hold on; plot(xlim,[0.5 0.5],'k--'); title('Accuracy');
set(gca,'XTickLabels',{'Congruent','Incongruent'});
hold on; scatter((1:2)+randn(sum(~nan_idx),1)./10,acc,50,'k','filled','MarkerFaceAlpha',0.3)

subplot(2,1,2);
barwitherr(squeeze(std(posterior))./sqrt(nSub),squeeze(mean(posterior)));
title('Classifier evidence'); legend('Imagined','Real');
set(gca,'XTickLabels',{'Congruent','Incongruent'});


%% Stats
if cfg.permutation

    pVals = nan(nTT,1);
    nBtstrp = 10000;

    for tt = 1:nTT
        bAcc = nan(nBtstrp,1);
        for b = 1:nBtstrp
            tmp = nan(size(acc,1),1);
            for sub = 1:size(acc,1)
                tmp(sub) = pAcc(sub,tt,randi(size(pAcc,3)));
            end
            bAcc(b) = mean(tmp);
        end

        pVals(tt) = sum(bAcc>squeeze(mean(acc(:,tt))))./nBtstrp;
    end

    % calculate the difference
    bAccDiff = nan(nBtstrp,1);
    for b = 1:nBtstrp
        tmp = nan(size(acc,1),1);
        for sub = 1:size(acc,1)
            rnum     = randi(size(pAcc,3));
            tmp(sub) = squeeze(pAcc(sub,1,rnum)-pAcc(sub,2,rnum));
        end
        bAccDiff(b) = mean(tmp);
    end
    pValDiff = sum(bAccDiff>mean(acc(:,1)-acc(:,2)))./nBtstrp;
    save(fullfile(group_dir,'ROI_decoding'),'pVals','pValDiff','-append');
end

%% Plot the trial distribution before and after
trl_distr(:,nan_idx,:,:) = [];
figure;
for b = 1:2
    dat = squeeze(trl_distr(b,:,:,:));
    M   = squeeze(mean(dat,1));
    SEM = squeeze(std(dat))./sqrt(size(dat,1));
    
    subplot(2,1,b); 
    barwitherr(SEM,M);
    set(gca,'XTickLabels',{'Congruent','Incongruent'});
end

%% Look at individual differences
cfg.plotIndividual = false; 
cfg.plotting = false;

[~,~,FA,H] = BehaviourAnalysis(cfg);
PP = (FA+H)./2;
shift = PP(:,1)-PP(:,2);
shift(nan_idx) = [];

split = shift<median(shift);
acc_split = cat(3,acc(split,:),acc(~split,:)); % low vs high

figure;
b = barwitherr(squeeze(std(acc_split))./sqrt(sum(split)),squeeze(mean(acc_split)));
hold on; plot(xlim,[0.5 0.5],'k--'); 
set(gca,'XTickLabels',{'Congruent','Incongruent'});
for m = 1:2
    hold on; scatter(b(m).XData+b(m).XOffset+randn(sum(~split),1)./20,squeeze(acc_split(:,:,m)),50,'k','filled','MarkerFaceAlpha',0.3)
end
ylim([0.15 0.8]);

% stats for individual differences
if cfg.permutation

    pValsSplit = nan(nTT,2);
    nBtstrp = 10000;

    pAccSplit = cat(4,pAcc(split,:,:),pAcc(~split,:,:));

    for tt = 1:nTT
        for m = 1:2
            bAcc = nan(nBtstrp,1);
            for b = 1:nBtstrp
                tmp = nan(size(acc_split,1),1);
                for sub = 1:size(acc_split,1)
                    tmp(sub) = pAccSplit(sub,tt,randi(size(pAccSplit,3)),m);
                end
                bAcc(b) = mean(tmp);
            end

            pValsSplit(tt,m) = sum(bAcc>squeeze(mean(acc_split(:,tt,m))))./nBtstrp;
        end
    end

    % calculate the difference    
    for m = 1:2
        bAccDiff = nan(nBtstrp,1);
        for b = 1:nBtstrp
            tmp = nan(size(acc_split,1),1);
            for sub = 1:size(acc_split,1)
                rnum     = randi(size(pAccSplit,3));
                tmp(sub) = squeeze(pAccSplit(sub,1,rnum,m)-pAccSplit(sub,2,rnum,m));
            end
            bAccDiff(b) = mean(tmp);
        end
        pValDiffSplit(m) = sum(bAccDiff>squeeze(mean(acc_split(:,1,m)-acc_split(:,2,m))))/nBtstrp;
    end
end
