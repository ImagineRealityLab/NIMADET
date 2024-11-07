function CheckMovement(cfg)
% function CheckMovement(cfg)

nSub = length(cfg.subjects);

% loop over subs
m_sub = nan(nSub,6,2); 
for sub = 1:nSub

    dr = dir(fullfile(cfg.root,cfg.subjects{sub},cfg.dir));
    dr = {dr.name};
    sess_dir = dr(contains(dr,'sess'));
    nSess = length(sess_dir);

    gcf = figure(1); m = nan(nSess,6,2);
    for s = 1:nSess
        nifti_dir = fullfile(cfg.root,cfg.subjects{sub},cfg.dir,sess_dir{s});
        movementFile = str2fullfile(nifti_dir,'rp*.txt');
        fid = fopen(movementFile);
        out = textscan(fid,'%f %f %f %f %f %f');
        %fclose(movementFile);
        movement = cell2mat(out); clear out
        m(s,:,1) = mean(movement,1);
        m(s,:,2) = max(abs(movement));

        subplot(2,nSess,s)
        plot(movement(:,1:3)); ylabel('mm'); xlabel('scans')
        legend('X','Y','Z');
        title(sprintf('%s session %d',cfg.subjects{sub},s))
        %ylim([-2 2])

        subplot(2,nSess,s+nSess)
        plot(movement(:,4:6)); ylabel('mm'); xlabel('scans')
        legend('pitch','roll','yaw');
        %ylim([-0.1 0.1])
    end
    set(gcf, 'Position', get(0, 'Screensize'));

    m_sub(sub,:,1) = squeeze(mean(m(:,:,1),1));
    m_sub(sub,:,2) = squeeze(max(m(:,:,2)));
    
    fprintf('Press any key to show next participant \n')
    pause;
end

figure(2);
for sub = 1:nSub
    subplot(5,6,sub)
    bar(squeeze(m_sub(sub,:,:)))
end




