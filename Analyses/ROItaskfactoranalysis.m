function act = ROItaskfactoranalysis(cfg)
% function ROItaskfactoranalysis(cfg)
%
% cfg.root = root;
% cfg.subjects = subjects(which_subjects);
% cfg.ROI      = fullfile(root,'Analyses','ROI_masks','amPFC_Simons.nii');
% cfg.data_dir = 'FirstLevel/MT_taskfactors_256_rmIC_2factors';
% cfg.cond{1}  = 'cong';
% cfg.cond{2}  = 'inco';

% load ROI mask
[~,mask] = read_nii(cfg.ROI);
nROIs    = length(unique(mask(:)))-1;

nSubs = length(cfg.subjects);
nCond = length(cfg.cond);

% get data
act   = nan(nSubs,nCond,nROIs);
for sub = 1:nSubs

    fprintf('Processing sub %d... \n',sub)
    FLdir = fullfile(cfg.root,'Results',cfg.subjects{sub},cfg.data_dir);

    load(fullfile(FLdir,'SPM.mat'),'SPM')

    nRuns = length(SPM.Sess);
    dat = nan(nCond,nRuns,nROIs);
    for cnd = 1:nCond
        for r = 1:nRuns
            beta_idx = find(contains(SPM.xX.name,cfg.cond{cnd}) & ...
            contains(SPM.xX.name,'bf(1)') & ...
            contains(SPM.xX.name,sprintf('Sn(%d)',r)));

            if ~isempty(beta_idx)
            [~,tmp] = read_nii(fullfile(FLdir,sprintf('beta_%04d.nii',beta_idx)));
            for roi = 1:nROIs
                dat(cnd,r,roi) = nanmean(tmp(mask(:)==roi));
            end
            end
        end
    end

    % mean centre per run
    if nRuns > 1
        zdat = nan(size(dat));
        for r = 1:nRuns
            if nROIs > 1
                zdat(:,r,:) = dat(:,r,:)-mean(dat(:,r,:),1);
            else
                zdat(:,r) = dat(:,r)-mean(dat(:,r));
            end
        end
    else
        zdat = dat;
    end

    % calculate mean activation
    %zdat = dat;
    act(sub,:,:) = squeeze(mean(zdat,2));
end

%% Save the data
if isfield(cfg,'outDir')
outDir  = fullfile(cfg.root,'Results',cfg.outDir);
if ~exist(outDir,'dir'); mkdir(outDir); end

save(fullfile(outDir,sprintf('%s.mat',cfg.ROI_name)),'act')
end
%% plot the results
if cfg.plot
figure(1); hold off;
if ~contains(cfg.cond{1},'x')
    act = act-mean(act,2); % mean centre per subject for non-pmod
end

if nROIs >1

    SEM = nan(nCond,nROIs); M = nan(nCond,nROIs);
    for r = 1:nROIs
        for c = 1:nCond
            dat = squeeze(act(:,c,r));
            nan_idx = isnan(dat);
            SEM(c,r) = std(dat(~nan_idx))./sqrt(sum(~nan_idx));
            M(c,r) = mean(dat(~nan_idx));
        end
    end
    b = barwitherr(SEM',M');

    if nCond == 4
        b(2).FaceColor = [0 0 0.5]; b(1).FaceColor = [0 0 1];
        b(4).FaceColor = [0.5 0 0]; b(3).FaceColor = [1 0 0];
    elseif nCond == 2
        b(1).FaceColor = [0 0 0.5];
        b(2).FaceColor = [0.5 0 0];
    end

    for c = 1:nCond
        dat = squeeze(act(:,c,:));
        hold on;
        %scatter(((1:nROIs)+b(c).XOffset)+(randn(nROIs,nSubs)./30)',dat,20,'k','filled','MarkerFaceAlpha',0.5)
    end
    set(gca,'XTickLabels',ROI_names);
else
    SEM = std(act)./sqrt(nSubs); M = mean(act);
    b   = barwitherr(SEM,M);
    hold on; scatter(b.XData+(randn(nCond,nSubs)./10)',act,40,'k','filled','MarkerFaceAlpha',0.3)
    set(gca,'XTickLabels',cfg.cond)
    title(cfg.ROI_name)
end
ylabel('% signal change')
set(gca,'Position',[0.1446 0.11 0.7606 0.8150])

%% Look at individual differences
cfg.plotIndividual = false; 
cfg.plotting = false;

[Cr,Dp,FA,H,V] = BehaviourAnalysis(cfg);
PP = (FA+H)./2;
shift = PP(:,1)-PP(:,2);
Dp = mean(Dp,2); C = mean(Cr,2); V = squeeze(mean(mean(V,2),3));
vals = [shift, Dp, C, V];
names = {'shift','Dp','C','V'};

if nROIs == 1 

    for v = 1:size(vals,2)
        
        % median split low vs high
        idx = vals(:,v)>median(vals(:,v)); 
        SEM = nan(size(act,2),2); M = nan(size(act,2),2);
        for i = 1:2
            SEM(:,i) = std(act(idx==(i-1),:))./sqrt(sum(idx==(i-1)));
            M(:,i)   = mean(act(idx==(i-1),:));
        end

        % plot all conditions
        figure(4);
        subplot(2,2,v);
        b = barwitherr(SEM,M);
        set(gca,'XTickLabels',{'cong pres','cong abs','inco pres','inco abs'});
        if v == 1; legend('low','high'); end
        act_split = cat(3,act(~idx,:),act(idx,:));
        for s = 1:2
            hold on; scatter(b(s).XData+b(s).XOffset+randn(sum(idx),1)./30,squeeze(act_split(:,:,s)),...
                40,'k','filled','MarkerFaceAlpha',0.3)
        end
        title(names{v});

        % plot cong versus inco
        cong = mean(act(:,1:2),2); inco = mean(act(:,3:4),2);
        dat = [cong, inco]; 
        figure(5)
        subplot(2,2,v)
        M = nan(2,2); SEM = nan(2,2);
        for i = 1:2
            M(i,:) = mean(dat(idx==(i-1),:));
            SEM(i,:) = std(dat(idx==(i-1),:))./sqrt(sum(idx==(i-1)));
        end        
        barwitherr(SEM',M')        
        set(gca,'XTick',1:2,'XTickLabels',{'Congruent','Incongruent'});
        if v == 1; legend('low','high'); end
        title(names{v});
    end

end


figure(7);
PP = (FA+H)./2;
diff = PP(:,1)-PP(:,2);
idx = diff>median(diff);

subplot(1,2,1);
cong = mean(act(:,1:2),2); inco = mean(act(:,3:4),2);
dat = [cong, inco]; 
M = nan(2,2); SEM = nan(2,2);
for i = 1:2
    M(i,:) = mean(dat(idx==(i-1),:));
    SEM(i,:) = std(dat(idx==(i-1),:))./sqrt(sum(idx==(i-1)));
end
hold off;
for i = 1:2
    errorbar([1 2],M(i,:),SEM(i,:),'Marker','.','MarkerSize',30); hold on
end
xlim([0.5 2.5]); legend('low','high')

subplot(1,2,2);
dat = squeeze(act(:,2)); % only cong FA's
M = nan(2,1); SEM = nan(2,1);
for i = 1:2
    M(i) = mean(dat(idx==(i-1))); 
    SEM(i) = std(dat(idx==(i-1)))./sqrt(sum(idx==(i-1)));
end
barwitherr(SEM,M); 
hold on; scatter((1:2)+randn(sum(idx),1)./10,[dat(~idx) dat(idx)],50,'k','filled','MarkerFaceAlpha',0.5)

cong = mean(act(:,[1 2]),2)-mean(act(:,[3 4]),2);
pres = mean(act(:,[1 3]),2)-mean(act(:,[2 4]),2);
int  = (act(:,1)-act(:,2))-(act(:,3)-act(:,4));
end