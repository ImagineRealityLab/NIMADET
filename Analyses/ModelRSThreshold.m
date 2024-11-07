function ModelRSThreshold(cfg)

subjects = cfg.subjects(which_subjects);
nSubs    = length(subjects);
nRuns    = 4;

outputDir = fullfile(cfg.root,'Results','GroupResults',cfg.outDir);
load(fullfile(outputDir,'BehaviouralModelFit'),'params');

%% Plot relationship between decision and RS for different conditions
RJ = cell(nSubs,1); RS = cell(nSubs,1);
Pres = cell(nSubs,1); Cond = cell(nSubs,1);
for sub = 1:nSubs
    
    % Get the behavioural data
    cfg.behDir = fullfile(fileparts(fileparts(cfg.root)),'Results',subjects{sub},'Regressors','Behaviour_matrix');   
    cfg.nRuns  = nRuns;
    [P,Vt,Vm,Dp,~,Cond{sub},Pres{sub}] = getBehaviour(cfg);

    % Use estimated parameters to generate predicted signals
    % get predicted signals
    [~, ~, RJ{sub}, ~, RS{sub}] = modelPredictions(P,Vt,Vm,Cond{sub},Pres{sub},Dp,params(sub,2,:),2,0);   

end

% plot data
figure; cs = ['b','r']; labs = {'Absent','Present'};
for p = 1:2
    subplot(2,2,p);
    for c = 1:2
        idx = Pres{2}==p-1 & Cond{2}==c;
        scatter(RS{2}(idx),RJ{2}(idx),40,cs(c),'filled','MarkerFaceAlpha',0.5); hold on
    end
    title(sprintf('Stimulus %s',labs{p}))
    xlabel('Reality signal'); ylabel('After thresholding')
end
legend('Congruent','Incongruent');

% calculate correlation between RS and RJ per cond
r_RSRJ = nan(nSubs,2,2);
for sub = 1:nSubs
    for p = 1:2
        for c = 1:2
            idx = Pres{sub}==(p-1) & Cond{sub}==c;
            r_RSRJ(sub,p,c) = corr(RS{sub}(idx),RJ{sub}(idx));
        end
    end
end

subplot(2,2,3:4)
b = barwitherr(squeeze(std(r_RSRJ))./sqrt(nSubs),squeeze(mean(r_RSRJ)));
ylim([0.95 1]); set(gca,'XTickLabels',{'Absent','Present'});
b(1).FaceColor = cs(1); b(2).FaceColor = cs(2);
