function VividnessRegressionDecodingSearchlight(cfg)
% trains and tests LDA classifiers with n-fold cross-validation on the task
% data

rng(1); % for reproduceability

lambda = 100000; % for ridge regularization 

[hdr,mask] = read_nii(cfg.func_mask);
mask = mask > 0.1;

[~,mind,cind] = searchlightIndices(mask,cfg.radius);
nSearchlights = length(mind);

nSub = length(cfg.subjects);

dbstop if warning

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

    X = data; nTrls = length(V); clear data;
    congruency = D==I;

    % Mean-center per run
    for r = 1:nRuns       
         X(run_idx==r,:) = X(run_idx==r,:)-mean(X(run_idx==r,:),1);

         % correct vividness for congruency and stimulus presence
         for c = 1:2
             for p = 1:2
                 idx = congruency==(c-1) & P==(p-1) & run_idx'==r;
                 V(idx) = V(idx)-mean(V(idx));
             end
         end              
    end

    %% Do the decoding per searchlight
    acc = zeros(hdr.dim);
    for s = 1:nSearchlights

        if s >= (nSearchlights/10) && mod(s,round((nSearchlights/10))) == 0
            fprintf('Progress: %d percent of searchlights \n',round((s/nSearchlights)*100))
        end

        x = X(:,mind{s});

        % remove nans
        nan_idx = any(isnan(x));
        x(:,nan_idx) = [];   

        % classify
        if ~isempty(x)            

            yhat     = nan(nTrls,1);                 

            % cross-validate based on runs
            RI       = unique(run_idx);
            for rn = 1:length(RI)
                run = RI(rn);

                train_idx = setdiff(1:nTrls,find(run_idx==run));
                trainx    = x(train_idx,:);
                testx     = x(run_idx==run,:);

                % ridge regression decoder
                I = eye(size(trainx, 2));     % Identity matrix (size of #predictors)
                beta_ridge = (trainx' * trainx + lambda * I) \ (trainx' * V(train_idx)');  % Ridge regression coefficients
                yhat(run_idx==run) = testx*beta_ridge;                            

            end

            % log the accuracy in centre voxel
            acc(cind{s}) = corr(yhat,V','Type','Spearman');
        end
    end

    % write results 
    outDir = fullfile(cfg.root,'Results',cfg.subjects{sub},cfg.outDir);
    if ~exist(outDir,'dir'); mkdir(outDir); end
    write_nii(hdr2,acc,fullfile(outDir,'accuracy.nii'));
    clear acc X; 
end

%% Group analysis
outDir = fullfile(cfg.data,'GroupResults',cfg.outDir);
if ~exist(outDir,'dir'); mkdir(outDir); end

acc = zeros([nSub hdr.dim]);
for sub = 1:nSub
    [~,tmp] = read_nii(fullfile(cfg.root,'Results',cfg.subjects{sub},cfg.outDir,'accuracy.nii'));
    acc(sub,:,:,:) = tmp; clear tmp
end

write_nii(hdr2,squeeze(mean(acc,1)),fullfile(outDir,'accuracy.nii'));

rpVals = zeros(hdr.dim);
[~,p] = ttest(acc(:,mask)); rpVals(mask) = 1-p;
write_nii(hdr2,rpVals,fullfile(outDir,'rpvals.nii'));


