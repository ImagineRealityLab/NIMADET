function VividnessRegressionDecodingROIs(cfg)
% trains and tests LDA classifiers with n-fold cross-validation on the task
% data

clearvars -except cfg

rng(1); % for reproduceability

cfg.nPerm = 25; % per participants
cfg.nBtstrp = 10000; % to create group null distribution

nROIs = length(cfg.ROIs);
for r = 1:nROIs
    [V,Y] = read_nii(fullfile(cfg.ROI_dir,[cfg.ROIs{r} '.nii']));
    if r == 1        
        atlas = zeros(V.dim);
    end
    atlas(Y(:)>0) = r; clear V Y
end

%[V,atlas] = read_nii('D:\NIMADET\Data\rvisual_ROIs_Kastner.nii');

mask      = atlas>0; atlas = atlas(mask);
nSub = length(cfg.subjects);

acc = nan(nSub,nROIs);
pAcc = nan(nSub,nROIs,cfg.nPerm);
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
        data = cat(1,data,betas); clear betas %
    end

    X = data; nTrls = length(V);

    % Mean-center per run
    for r = 1:nRuns       
         X(run_idx==r,:) = X(run_idx==r,:)-mean(X(run_idx==r,:),1);
         V(run_idx==r) = V(run_idx==r)-mean(V(run_idx==r));
    end

    %% Do the decoding per ROI per trial type
    for r = 1:nROIs

        fprintf('Decoding ROI %d out of %d \n',r,nROIs)

        x = X(:,atlas==r);

        % remove nans
        nan_idx = any(isnan(x));
        x(:,nan_idx) = [];   

        % classify
        if ~isempty(x)            

            yhat     = nan(nTrls,1); post = nan(nTrls,length(unique(V)));
            yhatP    = nan(nTrls,cfg.nPerm);                   

            % cross-validate based on runs
            RI       = unique(run_idx);
            for rn = 1:length(RI)
                run = RI(rn);

                train_idx = setdiff(1:nTrls,find(run_idx==run));

                % PCA
                warning('off')
                [coeff,score,~,~,exp] = pca(x(train_idx,:));
                idx = find(cumsum(exp)>80);                
                trainx = score(:,1:idx(1));
                testx  = x(run_idx==run,:)*coeff; 
                testx = testx(:,1:idx(1));

                % linear decoder
                weights = regress(V(train_idx)',trainx);
                yhat(run_idx==run) = testx*weights;       
                                
                % permutation
                for p = 1:cfg.nPerm
                    permY   = V(train_idx);
                    permY = permY(randperm(length(permY)));
                    weights = regress(permY',trainx);
                    yhatP(run_idx==run,p) = testx*weights;
                end
                warning('on')
            end

            acc(sub,r) = corr(yhat,V','Type','Spearman');
            for p = 1:cfg.nPerm   
                pAcc(sub,r,p) = corr(yhatP(:,p),V','Type','Spearman');
            end
        end
    end
end

%% Group analysis
outDir = fullfile(cfg.data,'GroupResults',cfg.outDir);
if ~exist(outDir,'dir'); mkdir(outDir); end

% remove nans
nan_idx = any(isnan(acc'));
acc(nan_idx,:) = []; pAcc(nan_idx,:,:) = []; nSub = sum(~nan_idx);

pVals = nan(nROIs,1);
for r = 1:nROIs
    bAcc = nan(cfg.nBtstrp,1);
    for b = 1:cfg.nBtstrp
        tmp = nan(nSub,1);
        for sub = 1:nSub
            tmp(sub) = pAcc(sub,r,randi(cfg.nPerm));
        end
        bAcc(b) = mean(tmp); clear tmp;
    end
    pVals(r) = sum(bAcc>squeeze(mean(acc(:,r),1)))/cfg.nBtstrp;
end

save(fullfile(outDir,'ROIdecoding.mat'),'acc','pAcc','pVals','cfg');


%% Plot results
figure;
b = barwitherr(std(acc)./sqrt(nSub),mean(acc));
set(gca,'XTickLabels',cfg.ROIs); 
