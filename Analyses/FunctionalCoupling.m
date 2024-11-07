function FunctionalCoupling(cfg)

nROIs = length(cfg.ROI);
nSubs = length(cfg.subjects);

% get masks
ROIs = cell(nROIs,1);
for r = 1:nROIs
    [~,ROIs{r}] = read_nii(cfg.ROI{r});
end

r_vals = nan(nSubs,nROIs,nROIs,2,2); % per congruency per presence

% Loop over subjects
cP = 1; cI = 2; cPr = 3; cR = 4; cV = 5; cIC = 6;
for sub = 1:nSubs

    fprintf('Processing subject %d out of %d \n',sub,nSubs)

    % behavioural data
    beh_dir = fullfile(cfg.root,'Results',cfg.subjects{sub},'Regressors','Behaviour_matrix');
    beh_files = str2fullfile(beh_dir,'run*.mat');

    % get the brain data
    beh_data = [];
    for b = 1:length(beh_files)
        load(beh_files{b},'B');
        beh_data = cat(1,beh_data,B); clear B
    end
    beh_data(beh_data(:,cIC)==0,:) = []; % remove incorrect

    % brain data
    dataDir = fullfile(cfg.root,'Results',cfg.subjects{sub},cfg.data);
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
        load(fullfile(cfg.root,'Data',cfg.subjects{sub},...
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

        nTrls = length(beta_names);
        act   = nan(nTrls,nROIs);
        for t = 1:nTrls
            [~,tmp] = read_nii(nifti_names{t});
            for r = 1:nROIs
                act(t,r) = nanmean(tmp(ROIs{r}(:)>0));
            end

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
        data = cat(1,data,act-mean(act,1)); % mean center per run
    end

    % define congruency
    congruency = D ~= I; % code incongruent, so congruent is first (0)

    % look at connectivity for congruency
    for c = 1:2
        for pres = 1:2
            idx = (congruency==(c-1)) & P==(pres-1);
            if sum(idx)>2
                r_vals(sub,:,:,c,pres) = corr(data(idx,:),data(idx,:));
            end
        end
    end

end

%% Save 
outputDir = fullfile(cfg.root,'Results',cfg.outDir);
if ~exist(outputDir,'dir'); mkdir(outputDir); end

save(fullfile(outputDir,'Results'),"r_vals",'cfg');

%% Plot connectivity per condition per ROI
[~,ROI_names] = fileparts(cfg.ROI);

% FG to rest
figure;
dat = squeeze(r_vals(:,1,2:end,:,:));
ROI_names = ROI_names(2:end);
for r = 1:(nROIs-1)
    subplot(ceil((nROIs-1)/2),2,r);

    tmp = squeeze(dat(:,r,:,:));

    M   = squeeze(mean(tmp,1)); 
    SEM = squeeze(std(tmp))./sqrt(nSubs);   

    b = barwitherr(SEM',M'); 
    for c = 1:2
        hold on; scatter(b(c).XData+b(c).XOffset+randn(nSubs,1)./20,...
            squeeze(tmp(:,c,:)),30,'k','filled','MarkerFaceAlpha',0.3)
    end
    b(1).FaceColor = 'b'; b(2).FaceColor = 'r';
    set(gca,'XTickLabels',{'Absent','Present'})
    title(ROI_names(r))
end
