function [RJ,Vt,Vm,Dp,Cr,Cond,Pres] = getBehaviour(cfg)

B = [];
for r = 1:cfg.nRuns
    dat = load(fullfile(cfg.behDir,sprintf('run_%d.mat',r)),'B');
    B   = cat(1,B,dat.B); clear dat
end

% remove incorrect blocks and any nans
B(B(:,6)==0,:) = [];
nan_idx = any(isnan(B'));
B(nan_idx,:) = [];

% reorganize to cond (1 = cong, 2 = inco) x present x RJ x vividness 
Cond       = double(B(:,1)~=B(:,2))+1;
Pres       = B(:,3);
RJ         = B(:,4);
Vt         = B(:,5);

% get mean values
load(fullfile(fileparts(fileparts(cfg.behDir)),'behResults'));
Vm = squeeze(mean(Vt));
Dp = squeeze(mean(results.D));
Cr = squeeze(mean(results.C));

