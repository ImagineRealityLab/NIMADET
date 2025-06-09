function RSPredictedAct(cfg)
% checks whether the MV pattern underlying vividness representations also
% follows the other effects of a RS as predicted by model simulations: main
% effect of congruency and stimulus presence and vividness still increasing
% RS even within RJ's

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
nSub = length(cfg.subjects);
RS_congstim = nan(nSub,nROIs,2,2);
RS_vivresp  = nan(nSub,nROIs,4); 
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
    
    congruency = D~=I; % code for congruent first  

    X = data; nTrls = length(V); clear data;

    % Mean-center per run
    for r = 1:nRuns       
         X(run_idx==r,:) = X(run_idx==r,:)-mean(X(run_idx==r,:),1);
         V(run_idx==r) = V(run_idx==r)-mean(V(run_idx==r));
    end

    %% Do the decoding per searchlight
    for r = 1:nROIs

        x = X(:,atlas==r);

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

            % log the predicted RS for each condition            
            for c = 1:2
                for s = 1:2
                    idx = congruency==(c-1) & P==(s-1); % absent first 
                    RS_congstim(sub,r,c,s) = mean(yhat(idx));
                end
            end

            % log predicted RS for each response type
            RS_vivresp(sub,r,1) = mean(yhat(R==0 & V<(median(V(R==0)))));
            RS_vivresp(sub,r,2) = mean(yhat(R==0 & V>(median(V(R==0)))));
            RS_vivresp(sub,r,3) = mean(yhat(R==1 & V<(median(V(R==1)))));
            RS_vivresp(sub,r,4) = mean(yhat(R==1 & V>(median(V(R==1)))));
           
        end
    end
end

%% Plot the results
figure;
for r = 1:nROIs
    
    tmp1 = squeeze(RS_congstim(:,r,:,:));
    tmp2 = squeeze(RS_vivresp(:,r,:));

    subplot(2,nROIs,r) 
    M = squeeze(mean(tmp1));
    SEM = squeeze(std(tmp1))./sqrt(nSub);
    b = barwitherr(SEM',M'); title(cfg.ROIs{r});
    for c = 1:2
        hold on; scatter(b(c).XData+b(c).XOffset+randn(nSub,1)./30,squeeze(tmp1(:,c,:)),30,'k','filled','MarkerFaceAlpha',0.3)
    end
    set(gca,'XTickLabels',{'Absent','Present'});   

    subplot(2,nROIs,r+nROIs)
    M = squeeze(mean(tmp2));
    SEM = squeeze(std(tmp2))./sqrt(nSub);
    barwitherr(SEM,M); 
    hold on; scatter((1:4)+randn(nSub,1)./10,tmp2,30,'k','filled','MarkerFaceAlpha',0.3)

end
