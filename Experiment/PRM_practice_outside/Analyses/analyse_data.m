function analyse_data(cfg)

nsubs = length(cfg.files);

% loop over data-sets
for s = 1:nsubs
    
    % load the data
    dataFile = str2fullfile(fullfile(cfg.dir,cfg.files{s}),'PMT*');
    load(dataFile{1},'R','C','blocks','miniblocks','trials')
    
    % determine congruency
    nMB = length(miniblocks);
    ima = reshape(repmat(blocks',cfg.nMB,1),nMB,1);
    det = miniblocks;
    cong = ima==det;
    
    % get acc per mini-block
    acc = nan(nMB,1);
    for mb = 1:nMB
        acc(mb) = mean(squeeze(R(mb,:,3))==squeeze(trials(mb,:)));
    end
    
    % plot
    if cfg.plot
        figure(1);
        subplot(2,1,1);
        bar([mean(acc(cong)) mean(acc(~cong))]);
        set(gca,'XTickLabel',{'Congruent','Incongruent'})
        ylabel('Accuracy');

        subplot(2,1,2);
        viv = mean(squeeze(R(:,:,1)),2);
        bar([mean(viv(cong)) mean(viv(~cong))]);
        set(gca,'XTickLabel',{'Congruent','Incongruent'})
        ylabel('Vividness');
    end
    
    fprintf('\t Accuracy for %s is: \n',cfg.files{s})
    acc'
    fprintf('\t Mean vividness is:')
    mean(squeeze(R(:,:,1)),2)'

end
