% executes all sub-scripts in the correct order and collects intermediate
% results
clear; plotting = false;
addpath('Utilities');

% settings
orientations = [135 75 15]; % rotAngles of gratings to use
orientations = orientations(randperm(length(orientations)));
orientations = orientations(1:2); % pick two randomly
vivDet       = 2; % ALWAYS DETECTION FIRST
subID        = input('Please enter participant name: ','s'); % get the subject ID
outDir       = fullfile('Results',subID); 
if ~exist(outDir,'dir'); mkdir(outDir); else; warning('Directory already exists!'); end

save(fullfile(outDir,sprintf('%s_settings',subID)),'orientations','vivDet')

fprintf('\t THE ORIENTATIONS ARE %d AND %d \n \t VIVDET IS %d \n',orientations,vivDet);

% =========================================================================
% Practice detection - outside scanner
% =========================================================================
% just a few trials per orientation to get an idea of what gratings in 
% noise look like 

[PD,PA,V] = practiceDetection(subID,orientations); % also gives a rough 
% initial estimate of V (visibility)

% =========================================================================
% Initial staircase - outside scanner
% =========================================================================
% Provide an initial estimate of the contrast, to be fine tuned in the
% scanner later

[SC_V1,SC_acc1] = outsideStaircase(subID,orientations(1),'A',V,'test'); % grating 1
save(fullfile(outDir,sprintf('%s_settings',subID)),'SC_V1','-append')

[SC_V2,SC_acc2] = outsideStaircase(subID,orientations(2),'B',V,'test'); % grating 2
save(fullfile(outDir,sprintf('%s_settings',subID)),'SC_V2','-append')

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

[~,b1] = min(abs(SC_acc1-0.7));
[~,b2] = min(abs(SC_acc2-0.7));
V   = mean([SC_V1(b1) SC_V2(b2)]); % take mean staircased SAME FOR BOTH STIMULI!?

% =========================================================================
% Practice imagery - outside scanner
% =========================================================================
% includes a pretty difficult discimrination task, practice is done once
% participant has 10 correct trials ina row

[PI_V,PI_P] = practiceImagery(subID,orientations);

% =========================================================================
% Practice main task - outside scanner
% =========================================================================
% includes all elements of the main task, working at a rough estimate of
% the visibility                                                 % visibility to practice task 
PMT = practiceMainTask(subID,orientations,vivDet,V);
sca