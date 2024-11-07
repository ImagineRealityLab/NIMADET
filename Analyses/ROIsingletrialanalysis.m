function act = ROIsingletrialanalysis(cfg)
% function ROIsingletrialanalysis(cfg)
%
% cfg = [];
% cfh.root = root;
% cfg.subjects = subjects(which_subjects);
% cfg.ROI      = fullfile(root,'Analyses','ROI_masks','amPFC_Simons.nii');
% cfg.data_dir = 'FirstLevel/MT_perTrialPoldrack';
% cfg.cond{1}  = 'P==1';
% cfg.cond{2}  = 'P==0';

% load ROI mask
[~,mask] = read_nii(cfg.ROI);
nROIs    = length(unique(mask(:)))-1;
ROI_names = {'midOcc'};

% exclude s33
% if any(contains(cfg.subjects,'S33'))
%     tmp = find(contains(cfg.subjects,'S33'));
%     cfg.subjects(tmp)=[];
% end

% get data
nSubs = length(cfg.subjects);
nRuns = 4;
nCnds = length(cfg.cond);

act   = nan(nSubs,nCnds,nROIs);
for sub = 1:nSubs

    % get the trial indices per run per condition
    trl_idx = cell(nRuns,nCnds);
    for r = 1:nRuns
        load(fullfile(cfg.root,'Results',cfg.subjects{sub},...
            cfg.beh_dir,sprintf('run_%d.mat',r)),'B')
        nTrls = length(B);

        P = B(:,3); R = B(:,4); V = B(:,5);
        Pori = B(:,1); Iori = B(:,2);
        cong = B(:,1)==B(:,2);
        for cnd = 1:nCnds
            trl_idx{r,cnd} = eval(cfg.cond{cnd});
        end
        clear B P R cong
    end

    % load the data
    fprintf('Loading data for sub %d... \n',sub)
    dat = nan(nRuns,nTrls,nROIs);
    for r = 1:nRuns
        for t = 1:nTrls
            tmp  = str2fullfile(fullfile(cfg.root,'Results',cfg.subjects{sub},...
                cfg.data_dir,sprintf('RUN_%02d',r)),sprintf('trl_%d_*.nii',t));
            if ~isempty(tmp)
                [~,tmp] = read_nii(tmp);
                for roi = 1:nROIs
                    dat(r,t,roi) = nanmean(tmp(mask(:)==roi));
                end
            end
            clear tmp
        end
    end

    % mean center per run
    zdat = nan(size(dat));
    for r = 1:nRuns
        if nROIs ==1
            zdat(r,:) = dat(r,:)-mean(dat(r,:),2);
        else
            zdat(r,:,:) = dat(r,:,:)-mean(dat(r,:,:),2);
        end
    end

    % get mean per condition
    for c = 1:nCnds
        if nROIs == 1
            dat = [];
            for r = 1:nRuns
                dat = cat(1,dat,squeeze(zdat(r,trl_idx{r,c}))');
            end
            nan_idx = isnan(any(dat,2));
            act(sub,c) = mean(dat(~nan_idx,:));
        else
            for roi = 1:nROIs
                dat = [];
                for r = 1:nRuns
                    dat = cat(1,dat,squeeze(zdat(r,trl_idx{r,c},roi))');
                end
                nan_idx = isnan(any(dat,2));
                act(sub,c,roi) = mean(dat(~nan_idx,:));
            end
        end
    end

end

%% Remove nans
nan_idx = any(isnan(act'));
act(nan_idx,:) = []; 
nSubs = length(act);

%% Plot the activation
figure;
SEM = nan(nCnds,nROIs); M = nan(nCnds,nROIs);
for r = 1:nROIs
    for c = 1:nCnds
        dat = squeeze(act(:,c,r));
        nan_idx = isnan(dat);
        SEM(c,r) = std(dat(~nan_idx))./sqrt(sum(~nan_idx));
        M(c,r) = mean(dat(~nan_idx));
    end
end
barwitherr(SEM',M');
legend(cfg.cond)
set(gca,'XTickLabel',ROI_names)
hold on; scatter((1:4)+randn(nSubs,1)./10,act,30,'k','filled','MarkerFaceAlpha',0.5)

%% Do linear regression
Ycong = cat(1,act(:,1),act(:,2),act(:,3),act(:,4));
Xcong = [ones(length(Ycong),1), reshape(repmat([1:4],nSubs,1),nSubs*4,1)];

[betas,CI] = regress(Ycong,Xcong);

betas = nan(nSubs,2,2);
X     = [ones(4,1) (1:4)'];
for sub = 1:nSubs
    betas(sub,1,:) = regress(act(sub,1:4)',X);
    betas(sub,2,:) = regress(act(sub,5:8)',X);
end

SEM = squeeze(std(betas))./sqrt(nSubs);
M   = squeeze(mean(betas));

figure;
barwitherr(SEM',M')
