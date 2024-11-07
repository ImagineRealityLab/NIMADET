% executes all sub-scripts in the correct order and collects intermediate
% results
clear; 
addpath('Utilities');

% settings
subID        = input('Please enter participant name: ','s'); % get the subject ID
outDir       = fullfile('Results',subID); 
vivDet       = 2; % ALWAYS DETECTION FIRST

% load settings
load(fullfile(outDir,sprintf('%s_settings.mat',subID)),...
    'orientations','SC_V1','SC_V2');

% =========================================================================
% Staircase visibil71191ity - in scanner
% =========================================================================
 
V2_start = SC_V2(end);
[SC_V2,SC_acc2] = insideStaircase(subID,orientations(2),'B',V2_start,'mri'); % grating 2
pause;

% plot
subplot(2,1,1);
plot(SC_V1,'-o'); hold on; plot(SC_V2,'-o');
legend('Orientation 1','Orientation 2');
ylabel('Visibility')

subplot(2,1,2);
plot(SC_acc1,'-o'); hold on; plot(SC_acc2,'-o');
hold on; plot(xlim,[0.7 0.7],'k--')
legend('Orientation 1','Orientation 2');
ylabel('Accuracy'); 

% determine staircased values
V1 = SC_V1(end);
V2 = SC_V2(end);
[~,b1] = min(abs(SC_acc1-0.7));
[~,b2] = min(abs(SC_acc2-0.7));
V   = mean([SC_V1(b1) SC_V2(b2)]); % use same throughout

% =========================================================================
% Main task - in scanner
% =========================================================================
responseMappings = [1 1 2 2]; responseMappings = responseMappings(randperm(4));
save(fullfile(outDir,sprintf('RMs_%s',subID)),'responseMappings');
KbQueueRelease;


[R,C,T] = mainTask(subID,orientations,[V V],vivDet,'mri',1,responseMappings(1));
save(fullfile(outDir,sprintf('main_%s_run_%d',subID,1)),'R','C','T');
pause;

[R,C,T] = mainTask(subID,orientations,[V1 V2],vivDet,'mri',2,responseMappings(2));
save(fullfile(outDir,sprintf('main_%s_run_%d',subID,2)),'R','C','T');
pause;

[R,C,T] = mainTask(subID,orientations,[V1 V2],vivDet,'mri',3,responseMappings(3));
save(fullfile(outDir,sprintf('main_%s_run_%d',subID,3)),'R','C','T');
pause;

[R,C,T] = mainTask(subID,orientations,[V1 V2],vivDet,'mri',4,responseMappings(4));
save(fullfile(outDir,sprintf('main_%s_run_%d',subID,4)),'R','C','T');
pause;

% =========================================================================
% Localizer in the scanner
% =========================================================================
TR = 1.45;
orientations = [135 75 15]; % localizer of all orientations

runLocalizer(subID,orientations,TR,'mri',1)
pause;sca

runLocalizer(subID,orientations,TR,'mri',2)
pause;
