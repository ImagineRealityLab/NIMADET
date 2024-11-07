function [ w, rect ] = setWindow( debug )
%open psychtoolbox and set up screen for experiment
resolution = [];
if debug    
    resolution = [750 50 1250 550];
    %PsychDebugWindowConfiguration()     
end


% Run tests in experiment, skip during debugging
Screen('Preference','SkipSyncTests', 1)
screens = Screen('Screens');
screenNumber = max(screens);
doublebuffer = 1;

% Open screen with grey background
[w, rect] = Screen('OpenWindow', screenNumber,...
    [255/2,255/2,255/2],resolution, 32, doublebuffer+1);

% Set useful paramaters
KbName('UnifyKeyNames');
AssertOpenGL;
PsychVideoDelayLoop('SetAbortKeys', KbName('Escape'));
HideCursor();
Priority(MaxPriority(w));
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%The fMRI button box does not work well with KbCheck. I use KbQueue
%instead here, to get precise timings and be sensitive to all presses.
% KbQueueCreate;
% KbQueueStart;

end

