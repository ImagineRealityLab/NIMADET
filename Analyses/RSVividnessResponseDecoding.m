function RSVividnessResponseDecoding(cfg)
% checks whether the MV pattern underlying vividness representations can
% predict reality judgements, seperately for congruent and incongruent
% conditions (to check if the pattern is sensitive to stimulus content). To
% ensure differences between congruent and incongruent conditions are not
% due to differences in behavioural correlations (more high vividness
% ratings in real judgements in congruent vs incongruent conditions), the
% behavioural distributions are first equalized through downsampling

rng(1); % for reproduceability
lambda = 100000; % for ridge regularization 

% get ROIs
group_dir = fullfile(cfg.root,'GroupResults',cfg.dir);
nROIs     = length(cfg.ROIs);
for r = 1:nROIs
    [~,ROI] = read_nii(fullfile(group_dir,[cfg.ROIs{r} '.nii']));
    if r ==1
        atlas = ROI; 
    else
        atlas = atlas+(ROI.*r); 
    end
    clear ROI
end
mask = atlas>0; atlas = atlas(mask);

% loop over subjects
nEq  = 100; % repeat equalization procedure n times
nSub = length(cfg.subjects);
acc     = nan(nSub,nROIs,2,2,nEq); % per congruency, vividness and RJ
for sub = 1:nSub

    fprintf('Processing subject %s \n',  cfg.subjects{sub})

    % task data
    dataDir = fullfile(cfg.data,cfg.subjects{sub},'FirstLevel',['MT' cfg.dataDir]);
    nRuns   = length(dir(dataDir))-2;
    D = []; I = []; P = []; R = []; V = []; data = []; run_idx = [];
    for run = 1:nRuns

        fprintf('Reading in task betas for run %d... \n',run)

        if contains(cfg.dataDir,'Poldrack')
            beta_niftis = str2fullfile(fullfile(dataDir,sprintf('RUN_%02d',run)),'trl_*.nii');
            nifti_names = cell(length(beta_niftis),1); % re-sort
            beta_names  = cell(length(nifti_names),1);
            for n = 1:length(beta_niftis)
                nifti_names{n} = str2fullfile(fullfile(dataDir,sprintf('RUN_%02d',run)),sprintf('trl_%d_*.nii',n));
                [~,beta_names{n}] = fileparts(nifti_names{n});
            end
        else
            load(fullfile(dataDir,sprintf('RUN_%02d',run),'SPM.mat'),'SPM');
            bidx = find(contains(SPM.xX.name,'bf(1)') & contains(SPM.xX.name,'Ima'));
            beta_niftis = str2fullfile(fullfile(dataDir,sprintf('RUN_%02d',run)),'beta_*.nii');
            beta_names = SPM.xX.name(bidx);
            nifti_names = beta_niftis(bidx);
        end

        nTrls = length(beta_names);

        % remove blocks failed imagery check
        if cfg.ima_check
            load(fullfile(fileparts(cfg.data),'Data',cfg.subjects{sub},...
                'Behaviour',sprintf('MT_%s_run%d.mat',cfg.subjects{sub},run)),'C')
            if any(C==0)
                ridx = [];
                for c = 1:length(C)
                    idx = (c-1)*(nTrls/length(C))+1:c*(nTrls/length(C));
                    if C(c)==0
                        ridx = [ridx,idx];
                    end
                end
                beta_names(ridx) = [];
                nifti_names(ridx) = [];
            end
        end

        nTrls = length(beta_names);
        betas = nan(nTrls,sum(mask(:)));
        for t = 1:nTrls
            [hdr2,Y] = read_nii(nifti_names{t});
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
        data = cat(1,data,betas); clear betas %
    end
    
    congruency = D==I;  

    for eq = 1:nEq

        % Equalize per run
        X = []; V_eq = []; R_eq = []; run_idx_eq = []; congruency_eq = [];
        RI = unique(run_idx); nRuns = length(RI);
        for rn = 1:nRuns
            run = RI(rn);

            % Chart distributions
            v_vals = unique(V); nV = length(v_vals);
            old_distr_N = nan(2,2,nV); % per c, per rj, per v
            old_distr_idx = cell(2,2,nV);
            for c = 1:2
                for r = 1:2
                    for v = 1:nV
                        idx = run_idx'==run & congruency==(c-1) & R==(r-1) & V==v;
                        old_distr_idx{c,r,v} = find(idx);
                        old_distr_N(c,r,v)   = sum(idx);
                    end
                end
            end

            % Equalize vividness and RJs per congruency
            new_distr_idx = cell(2,2,nV);
            for r = 1:2
                for v = 1:nV
                    cong = old_distr_N(1,r,v);
                    inco = old_distr_N(2,r,v);
                    if cong==0 || inco==0 % if either is zero
                        new_distr_idx{1,r,v} = []; % empty
                        new_distr_idx{2,r,v} = [];
                    elseif cong > inco % downsample cong
                        rnd_idx = randperm(cong);
                        new_distr_idx{1,r,v} = old_distr_idx{1,r,v}(rnd_idx(1:inco));
                        new_distr_idx{2,r,v} = old_distr_idx{2,r,v};
                    elseif inco > cong % downsample inco
                        rnd_idx = randperm(inco);
                        new_distr_idx{2,r,v} = old_distr_idx{2,r,v}(rnd_idx(1:cong));
                        new_distr_idx{1,r,v} = old_distr_idx{1,r,v};
                    else % already equal
                        new_distr_idx{1,r,v} = old_distr_idx{1,r,v};
                        new_distr_idx{2,r,v} = old_distr_idx{2,r,v};
                    end
                end
            end

            % downsample
            equalized_idx = cell2mat(new_distr_idx(:)');
            X = cat(1,X,data(equalized_idx,:));
            V_eq = cat(1,V_eq,V(equalized_idx)');
            R_eq = cat(1,R_eq,R(equalized_idx)');
            run_idx_eq = cat(1,run_idx_eq,run_idx(equalized_idx));
            congruency_eq = cat(1,congruency_eq,congruency(equalized_idx)');
        end

        % Mean-center per run
        for r = 1:nRuns
            X(run_idx_eq==r,:) = X(run_idx_eq==r,:)-mean(X(run_idx_eq==r,:),1);
            V_eq(run_idx_eq==r) = V_eq(run_idx_eq==r)-mean(V_eq(run_idx_eq==r));
            R_eq(run_idx_eq==r) = R_eq(run_idx_eq==r)-mean(R_eq(run_idx_eq==r));
        end

        %% Do the analyses per ROI
        for r = 1:nROIs

            % select ROI
            x = X(:,atlas==r);

            % remove nans
            nan_idx = any(isnan(x));
            x(:,nan_idx) = [];

            % Decode per congruency, cross-validate over runs
            for c = 1:2
                idx = congruency_eq==(c-1);
                x_cong = x(idx,:);
                r_cong = R_eq(idx);
                v_cong = V_eq(idx);
                rn_cong = run_idx_eq(idx);

                nTrls = sum(idx);
                yhat  = nan(nTrls,1);
                for rn = 1:nRuns
                    run = RI(rn);

                    test_idx = find(rn_cong==run);
                    train_idx = setdiff(1:nTrls,test_idx);

                    % train
                    beta_ridge = (x_cong(train_idx,:)' * x_cong(train_idx,:) + lambda * eye(size(x_cong(train_idx,:), 2))) \...
                        (x_cong(train_idx,:)' * v_cong(train_idx));  % Ridge regression coefficients

                    % get predicted
                    yhat(test_idx) = x_cong(test_idx,:)*beta_ridge;
                end

                % calculate accuracy
                acc(sub,r,c,1,eq) = corr(yhat,v_cong,'Type','Spearman');
                acc(sub,r,c,2,eq) = corr(yhat,r_cong,'Type','Spearman');
            end

        end
    end
end
mAcc = squeeze(mean(acc,5));

%% Plot the results
nan_idx = squeeze(any(any(any(isnan(mAcc),4),3),2));
mAcc(nan_idx,:,:,:) = []; nSub = sum(~nan_idx);

figure;
kinds = {'vividness','reality judgement'};
for a = 1:2
    subplot(2,1,a)
    M = squeeze(mean(mAcc(:,:,:,a)))/sqrt(nSub);
    SEM = squeeze(std(mAcc(:,:,:,a)));
    barwitherr(SEM,M);
    set(gca,'XTickLabels',cfg.ROIs);
    title(kinds{a}); ylabel('r');
end

%% Check variance over repetitions
r = 5;
tmp = squeeze(acc(~nan_idx,r,2,1,:)); % congruent vividness
figure;
for s = 1:20
    subplot(4,5,s);
    histogram(tmp(s,:));
end

