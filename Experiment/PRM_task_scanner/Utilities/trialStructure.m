function trials = trialStructure(miniblocks,nTrials)

nMB = size(miniblocks,1);
trials = nan(nMB,nTrials);

for mb = 1:nMB    
    
    trials(mb,:) = [zeros(1,nTrials/2) ones(1,nTrials/2)]';
    
    % shuffle per block
    trials(mb,:) = trials(mb,randperm(nTrials));     
    
end