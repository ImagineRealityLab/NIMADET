function VividnessDecodingROIs(cfg)
% trains and tests LDA classifiers with n-fold cross-validation on the task
% data

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
posterior = nan(nSub,nROIs,4,4);
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

    % Balance the vividness ratings per run
    unique_vals = unique(V);
    nVals = length(unique_vals);
    Y = []; X = []; rn_idx = []; trl_per_run = nan(nRuns,1);
    for rn = 1:nRuns
        V_run = V(run_idx==rn);
        X_run = data(run_idx==rn,:);

        if length(unique(V_run))==nVals
            count = histc(V_run,unique_vals);
            idx   = balance_trials(V_run,'downsample',min(count)); 
            trl_per_run(rn) = min(count);
            Y     = cat(1,Y,V_run(cell2mat(idx'))');
            X     = cat(1,X,X_run(cell2mat(idx'),:));
            rn_idx = cat(1,rn_idx,ones(length(cell2mat(idx')),1).*rn);
        end 
    end
    V = Y'; run_idx = rn_idx;
    nTrls = length(V);

%     % mean center per run
%     for r = 1:nRuns       
%         X(run_idx==r,:) = X(run_idx==r,:)-mean(X(run_idx==r,:),1);
%     end

    if length(unique(run_idx))>1 && ~(min(trl_per_run)==1 && length(unique(run_idx))<3) && length(unique_vals)==4


    %% Do the decoding per ROI per trial type
    for r = 1:nROIs

        fprintf('Decoding ROI %d out of %d \n',r,nROIs)

        x = X(:,atlas==r);

        % remove nans
        nan_idx = any(isnan(x));
        x(:,nan_idx) = [];

        % mean center over voxels
        % x = x-mean(x,2);

        % do PCA - ONLY ON TRAINING SET
%         warning('off');
%         [~,score,~,~,exp] = pca(x);
%         warning('on');
%         idx = find(cumsum(exp)>80);
%         x = score(:,1:idx(1));        

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
                weights = regress(V(train_idx)',[ones(length(train_idx),1) trainx]);
                yhat(run_idx==run) = [ones(sum(run_idx==run),1) testx]*weights;       
                %[yhat(run_idx==run),~,post(run_idx==run,:)] = classify(x(run_idx==run,:),x(train_idx,:),V(train_idx)','diagLinear');
                                
                % permutation
                for p = 1:cfg.nPerm
                    permY   = V(train_idx);
                    permY = permY(randperm(length(permY)));
                    weights = regress(permY',[ones(length(train_idx),1) trainx]);
                    yhatP(run_idx==run,p) = [ones(sum(run_idx==run),1) testx]*weights;
                    %yhatP(run_idx==run,p) = classify(x(run_idx==run,:),x(train_idx,:),permY','diagLinear');
                end
                warning('on')
            end

            %acc(sub,r) = mean(yhat==V');
            %v_vals = unique(V);
            %for v = 1:length(v_vals)
            %    idx = V'==v_vals(v);
            %    posterior(sub,r,v_vals(v),v_vals) = mean(post(idx,:));
            %end

            acc(sub,r) = corr(yhat,V','Type','Spearman');
            %[~,~,~,acc(sub,r)] = perfcurve(V,yhat,max(V));
            for p = 1:cfg.nPerm   
                %pAcc(sub,r,p) = mean(yhatP(:,p)==V');
                pAcc(sub,r,p) = corr(yhatP(:,p),V','Type','Spearman');
            end
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

posterior(nan_idx,:,:,:) = [];

%% Plot results
figure;
b = barwitherr(std(acc)./sqrt(nSub),mean(acc));
set(gca,'XTickLabels',cfg.ROIs); 
%ylim([0.2 0.35]); hold on; plot(xlim,[0.25 0.25],'k--')

figure;
for r = 1:nROIs
    subplot(3,4,r); 
    imagesc(squeeze(mean(posterior(:,r,:,:)))); caxis([0.1 0.4])
    title(cfg.ROIs{r})
end

roi = 2;
c_map = makeColorMaps('teals');
cs    = c_map(round(linspace(60,256,4)),:);
figure;
subplot(1,2,1); imagesc(squeeze(mean(posterior(:,roi,:,:)))); colorbar;
subplot(1,2,2);
M = squeeze(mean(posterior(:,roi,:,:)));
SEM = squeeze(std(posterior(:,roi,:,:)))./sqrt(size(posterior,1));
b = errorbar(M',SEM'); legend('1','2','3','4')
for v = 1:4
    b(v).Color = cs(v,:);
end

%% Test gradedness
% By testing whether the posterior probability for the neighbouring
% vividness rating is higher compared to the one next to that
neighbours(1,1) = 2;
neighbours(2,1) = 1;
neighbours(3,1) = 4;
neighbours(4,1) = 3;

neighbours(1,2) = 3;
neighbours(2,2) = 4;
neighbours(3,2) = 1;
neighbours(4,2) = 2;

roi = 2;
pp_N = nan(nSub,2); 
for sub = 1:nSub
   
    tmp = nan(nVals,2);
    for v = 1:nVals
        tmp(v,1) = squeeze(posterior(sub,roi,v,neighbours(v,1)));
        tmp(v,2) = squeeze(posterior(sub,roi,v,neighbours(v,2)));
    end
    pp_N(sub,:) = mean(tmp,1);
end

[h,p] = ttest(pp_N(:,1),pp_N(:,2));
