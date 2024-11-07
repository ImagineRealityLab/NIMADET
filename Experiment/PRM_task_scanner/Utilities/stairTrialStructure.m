function [ trialMatrix ] = stairTrialStructure(nTrials, nEvalTrials)
% Trial structure
% 50% present/50% absent

if nargin < 4
    nEvalTrials = 10;
end
nStairs = nTrials/nEvalTrials;

% Trial Matrix
trialMatrix = nan(nTrials, 1);

for s = 1:nStairs
    
    % staircase indices
    sIdx = (s-1)*nEvalTrials+1:s*nEvalTrials;
    tmp  = nan(nEvalTrials);
    
    % 50/50 absence presence
    tmp(1:nEvalTrials/2) = 0; 
    tmp(nEvalTrials/2+1:nEvalTrials) = 1;
    
    % shuffle within stairs
    trialMatrix(sIdx) = tmp(randperm(nEvalTrials));
end

end