function DecodeVividnessVersusRJ(cfg)
% function DecodeVividnessVersusRJ(cfg)
%
% Binarizes vividness ratings and then orthogonalizes vividness ratings and
% reality judgements to remove the behavioural correlation between them.
% Then applies LDA binary classification on both vividness and RJ to
% compare the accuracy between the two.

rng(1); % for reproduceability

ds = 25; % random downsample iterations 
cfg.nPerm = 25; % per participants
cfg.nBtstrp = 10000; % to create group null distribution

% get ROIs
group_dir = fullfile(cfg.root,'Results','GroupResults',cfg.dir);
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
nSub      = length(cfg.subjects);

acc   = nan(nSub,2,nROIs); % vividness and RJ
pAcc  = nan(nSub,2,nROIs,cfg.nPerm);

for sub = 1:nSub

    fprintf('Processing subject %s \n',  cfg.subjects{sub})

    % out dir
    outDir = fullfile(cfg.root,'Results',cfg.subjects{sub},cfg.outDir);
    if ~exist(outDir,'dir'); mkdir(outDir); end

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

    % balance trial numbers per run
    Vl = V>=median(V); % binarize
    ridx = unique(run_idx); nRuns = length(ridx);
    Y1 = []; Y2 = [];
    if cfg.balancetrialsRJViv % per RJ/Viv combination

        X = [];  Ridx = [];
        for r = 1:nRuns
            run = ridx(r);
            y1 = Vl(run_idx==run)'; % vividness
            y2 = R(run_idx==run)'; % RJ

            y3  = nan(length(y1),1);
            y3(y1==0&y2==0) = 1;
            y3(y1==0&y2==1) = 2;
            y3(y1==1&y2==0) = 3;
            y3(y1==1&y2==1) = 4;

            if length(unique(y3))==4 % only if we have all classes
                idx = balance_trials(y3,'downsample');
                Y1  = cat(1,Y1,y1(cell2mat(idx)));
                Y2  = cat(1,Y2,y2(cell2mat(idx)));
                X   = cat(1,X,data(cell2mat(idx),:));
                Ridx = cat(1,Ridx,ones(length(cell2mat(idx)),1)*run);
            end
        end
        run_idx = Ridx; clear Ridx;

    else % for smallest trial number, irrespective of V/RJ relationship

        X1 = []; X2 = [];
        Ridx1 = []; Ridx2 = [];
        for r = 1:nRuns
            run = ridx(r);
            y1 = Vl(run_idx==run)'; % vividness
            y2 = R(run_idx==run)'; % RJ

            nTrlsPerClass = [sum(y1==0) sum(y1==1) sum(y2==0) sum(y2==1)];
            nTrlsMin      = min(nTrlsPerClass);

            if ~any(nTrlsPerClass==0) % only if we have all classes
                % vividness
                idx1 = balance_trials(y1+1,'downsample',nTrlsMin);
                Y1   = cat(1,Y1,y1(cell2mat(idx1)));
                X1   = cat(1,X1,data(cell2mat(idx1),:));
                Ridx1 = cat(1,Ridx1,(ones(length(cell2mat(idx1)),1)*run));

                % RJ
                idx2 = balance_trials(y2+1,'downsample',nTrlsMin);
                Y2   = cat(1,Y2,y2(cell2mat(idx2)));
                X2   = cat(1,X2,data(cell2mat(idx2),:));
                Ridx2 = cat(1,Ridx2,(ones(length(cell2mat(idx2)),1)*run));
            end
        end
    end
    clear data; nTrls = length(Y2);

    if ~cfg.balancetrialsRJViv; X = X1; end
    if ~isempty(X)

        % mean center
        for r = 1:nRuns
            run = ridx(r);
            if cfg.balancetrialsRJViv
                X(run_idx==run,:) = X(run_idx==run,:)-nanmean(X(run_idx==run,:),2); % average over voxels
            else
                X1(Ridx1==run,:) = X1(Ridx1==run,:)-nanmean(X1(Ridx1==run,:),2); % average over voxels
                X2(Ridx2==run,:) = X2(Ridx2==run,:)-nanmean(X2(Ridx2==run,:),2); % average over voxels
            end
        end

        %% Do the decoding per ROI per judgement

        % loop over ROIs
        for r = 1:nROIs

            % get data for this ROI
            if cfg.balancetrialsRJViv
                x = X(:,atlas==r);
                nan_idx = any(isnan(x));
                x(:,nan_idx) = [];
            else
                x1 = X1(:,atlas==r);
                nan_idx = any(isnan(x1));
                x1(:,nan_idx) = [];
                x2 = X2(:,atlas==r);
                nan_idx = any(isnan(x2));
                x2(:,nan_idx) = [];
            end

            % cross-validate over runs
            if cfg.balancetrialsRJViv
                for rn = 1:nRuns
                    run = ridx(rn);

                    yhatR = nan(nTrls,1); yhatPR = nan(nTrls,cfg.nPerm);
                    yhatV = nan(nTrls,1); yhatPV = nan(nTrls,cfg.nPerm);
                    test_idx  = find(run_idx==run);
                    train_idx = setdiff(1:nTrls,test_idx);

                    % vividness
                    decoder         = train_LDA(cfg,Y1(train_idx),x(train_idx,:)');
                    yhatV(test_idx) = decode_LDA(cfg,decoder,x(test_idx,:)');

                    % reality judgement
                    decoder         = train_LDA(cfg,Y2(train_idx),x(train_idx,:)');
                    yhatR(test_idx) = decode_LDA(cfg,decoder,x(test_idx,:)');

                    % permutation
                    for p = 1:cfg.nPerm

                        % vividness
                        trainy = Y1(train_idx);
                        permY = trainy(randperm(length(trainy)));
                        decoder = train_LDA(cfg,permY,x(train_idx,:)');
                        yhatPV(test_idx,p) = decode_LDA(cfg,decoder,x(test_idx,:)');

                        % reality judgement
                        trainy = Y2(train_idx);
                        permY = trainy(randperm(length(trainy)));
                        decoder = train_LDA(cfg,permY,x(train_idx,:)');
                        yhatPR(test_idx,p) = decode_LDA(cfg,decoder,x(test_idx,:)');
                    end
                end

                % calculate accuracy
                acc(sub,1,r) = mean((yhatV>0)==Y1); pAcc(sub,1,r,:) = mean((yhatPV>0)==repmat(Y1,1,cfg.nPerm));
                acc(sub,2,r) = mean((yhatR>0)==Y2); pAcc(sub,2,r,:) = mean((yhatPR>0)==repmat(Y2,1,cfg.nPerm));

            else % balanced seperately
                for type = 1:2
                    yhat = nan(nTrls,1); yhatP = nan(nTrls,cfg.nPerm);
                    for rn = 1:nRuns
                        run = ridx(rn);
                        if type==1
                            test_idx  = find(Ridx1==run);
                            x         = x1;
                            y         = Y1;
                        elseif type==2
                            test_idx  = find(Ridx2==run);
                            x         = x2;
                            y         = Y2;
                        end

                        train_idx = setdiff(1:nTrls,test_idx);
                        decoder   = train_LDA(cfg,y(train_idx),x(train_idx,:)');
                        yhat(test_idx) = decode_LDA(cfg,decoder,x(test_idx,:)');

                        % permutation
                        for p = 1:cfg.nPerm
                            trainy = y(train_idx);
                            permY = trainy(randperm(length(trainy)));
                            decoder = train_LDA(cfg,permY,x(train_idx,:)');
                            yhatP(test_idx,p) = decode_LDA(cfg,decoder,x(test_idx,:)');
                        end
                    end
                    % calculate accuracy
                    acc(sub,type,r) = mean((yhat>0)==y); pAcc(sub,type,r,:) = mean((yhatP>0)==repmat(y,1,cfg.nPerm));
                end
            end
        end

    end
end

% save group-results
group_dir = fullfile(cfg.root,'Results','GroupResults',cfg.outDir);
if ~exist(group_dir,'dir'); mkdir(group_dir); end
save(fullfile(group_dir,'ROI_decodingNotBalancedRJViv'),'acc','pAcc');


% remove nans
nan_idx = squeeze(any(any(isnan(acc),2),3));
%fprintf('Removed %d nan subs \n',sum(nan_idx))
%load(fullfile(group_dir,'nan_idx.mat'),'nan_idx');
acc(nan_idx,:,:) = []; pAcc(nan_idx,:,:,:) = [];

%% Plot the results

figure;
b = barwitherr(squeeze(std(acc))'./sqrt(nSub),squeeze(mean(acc))');
ylim([0.4 0.6]); hold on; plot(xlim,[0.5 0.5],'k--'); title('Accuracy');
set(gca,'XTickLabels',cfg.ROIs)

for c = 1:2
    hold on;  scatter(b(c).XData+b(c).XOffset+randn(sum(~nan_idx),1)./20,squeeze(acc(:,c,:)),30,'k','filled','MarkerFaceAlpha',0.3)
end

%% Stats
if cfg.permutation

    pVals = nan(2,nROIs);
    nBtstrp = 10000;

    for tt = 1:2
        for r = 1:nROIs
            bAcc = nan(nBtstrp,1);
            for b = 1:nBtstrp
                tmp = nan(size(acc,1),1);
                for sub = 1:size(acc,1)
                    tmp(sub) = pAcc(sub,tt,r,randi(size(pAcc,3)));
                end
                bAcc(b) = mean(tmp);
            end

            pVals(tt,r) = sum(bAcc>squeeze(mean(acc(:,tt,r))))./nBtstrp;
        end
    end

    % difference
    pValDiff = nan(nROIs,1);
    for r = 1:nROIs
        bAcc = nan(nBtstrp,1);
        for b = 1:nBtstrp
            tmp = nan(size(acc,1),1);
            for sub = 1:size(acc,1)
                rndi = randi(size(pAcc,3));
                tmp(sub) = pAcc(sub,1,r,rndi)-pAcc(sub,2,r,rndi);
            end
            bAcc(b) = mean(tmp);
        end
        pValDiff(r) = sum(bAcc>(squeeze(mean(acc(:,1,r)-acc(:,2,r)))))./nBtstrp;
    end

end


