
%% outside staircase detection (SD) script
% Opens window -> defines trials -> makes gabors -> adds noise ->
% runs through trials.

function [V,acc] = outsideStaircase(subID,orientation,gratingName,V,environment)

% =========================================================================
% Setup
% =========================================================================

% threshold value
threshold = 0.7; % accuracy to staircase to

% output
output = fullfile(cd,'results',subID);
if ~exist(output,'dir'); mkdir(output); end
saveName = sprintf('SC_%d_%s.mat',orientation,subID);

[w, rect] = setWindow(0);
ShowCursor;

% Trial numbers
nEvalTrials = 6;                       % Determine detection accuracy for this many trials
nStairs = 8;
nTrials = nStairs*nEvalTrials;

% Get the trial structure
trials = stairTrialStructure(nTrials, nEvalTrials);

% Visibility settings
vis_scale = [0 logspace(log10(0.005),log10(0.2),299)]; % steps in log space
visibility = V;                         % Start visibility
[~,visibility] = min(abs(vis_scale-visibility));  % scale idx
lower_b ound = 0.6;                      % Go up around here
upper_bound = 0.8;                      % Go down around here
V = zeros(nStairs,1);                   % Track visibility
acc = zeros(nStairs,1);                 % Track detection acc

% responses
B = zeros(nTrials,2); 
trialResponse = 1;
trialRT       = 2;
if strcmp(environment,'mri')
    yesKey        = '1!';
    noKey         = '2@';
else
    yesKey        = 'h';
    noKey         = 'j';
end

% timing
fixTime       = 0.2;
mITI          = 2; % mean ITI - randomly sample from norm
sITI          = 1; % SD for sampling
ITIs          = normrnd(mITI,sITI,1,nTrials);

% =========================================================================
% Stimuli
% =========================================================================

% Makes the gabors to show for instruction
gaborPatch = make_stimulus(orientation,1); % full visibility
gaborTexture = Screen('MakeTexture',w,gaborPatch);

% Noise stimulus info
displayDuration = 2;                     % Duration of the stimulus in seconds
hz = Screen('NominalFrameRate', w);
ifi = 1/hz;                              % Refresh rate
nStepSize = 2;                           % 2 frames per step
nSteps = (displayDuration/ifi)/nStepSize;

%% Instructon screen
[xCenter, yCenter] = RectCenter(rect);
[x_pix, ~] = Screen('WindowSize', w);

% show instructions
text = ['We will now do a calibration on the visibility of the gratings in noise. \n ',...
    'Your task is simply to indicate whether the grating below is presented (Yes [RI] or No [RM]) \n ',...
    'The grating will be presented in 50% of the trials. Over the course of the block it \n ',...
    'will become harder to see the gratings. Try to focus as best as you can. \n ',...
    'Each block will take about 5 minutes. Good luck! \n \n ',...
    '[Press any key to start] \n '];

Screen('TextSize',w, 28);
DrawFormattedText(w, text, 'center', yCenter*0.75, [255 255 255]);

% show gratings
Xpos     = x_pix*(1/2);
baseRect = [0 0 size(gaborPatch,1) size(gaborPatch,1)];
allRect  = CenterRectOnPointd(baseRect,Xpos,yCenter*1.4);

Screen('DrawTextures', w, gaborTexture, [], allRect, [],[], 0.5);
DrawFormattedText(w, sprintf('Grating %s',gratingName), xCenter*0.94, yCenter*1.2, [255 255 255]);

vbl=Screen('Flip', w);
WaitSecs(1);
KbWait;
    
%% Trials start
stairCount = 0; stairStep = 1;
for iTrial = 1:nTrials
    
    % Fixation
    Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
        4, [0,0,0], [rect(3)/2, rect(4)/2], 1);
    vbl=Screen('Flip', w);
    WaitSecs(fixTime);
    
    % Grating noise frames
    if trials(iTrial) == 1 % present trial
        % schedule of visibility gradient (i.e how visible at each frame)
        % Increases till most visible at the end
        schedule = vis_scale(round(linspace(1,visibility,nSteps)));
    else % Pure noise trial
        % 0 for entire schedule
        schedule = zeros(1,nSteps);
    end
        
    % Make the texture for each frame by combining the gabor with noise.
    % Rotates the annulus mask to hide the rotated boundary box around the
    % grating.
    target = {};
    for i_frame = 1:nSteps
        idx = ((i_frame-1)*nStepSize)+1:(i_frame*nStepSize);
        tmp = Screen('MakeTexture',w, ...
            make_stimulus(orientation,schedule(i_frame)));
        for i = 1:length(idx); target{idx(i)} = tmp; end
    end
    
    % =========================================================================
    % Presentation
    % =========================================================================
    
    % Present stimulus
    tini = GetSecs;
    for i_frame = 1:length(target)
        
        while GetSecs-tini < ifi*i_frame
            
            Screen('DrawTextures',w,target{i_frame});
            
            Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
                4, [0,0,0], [rect(3)/2, rect(4)/2], 1); % fixation
            
            vbl=Screen('Flip', w);
        end
    end
    
    % Decision
    stairCount = stairCount + 1;
    text = 'Was there a grating on the screen? \n Yes [RI] or no [RM]';
    Screen('TextSize',w, 28);
    DrawFormattedText(w, text, 'center', 'center', 255);
    vbl = Screen('Flip', w);
    
    % Log response
    keyPressed = 0; % clear previous response
    while ~keyPressed
        
        [~, keyTime, keyCode] = KbCheck(-3);
        key = KbName(keyCode);
        
        if ~iscell(key) % only start a keypress if there is only one key being pressed
            if any(strcmp(key, {yesKey,noKey}))
                
                % fill in B
                B(iTrial,trialResponse) = strcmp(key,yesKey); % 1 yes 0 no
                B(iTrial,trialRT) = keyTime-vbl;
                
                keyPressed = true;
                
                elseif strcmp(key, 'ESCAPE')
                    Screen('TextSize',w, 28);
                    DrawFormattedText(w, 'Experiment was aborted!', 'center', 'center', [255 255 255]);
                    Screen('Flip',w);
                    WaitSecs(0.5);
                    ShowCursor;
                    disp(' ');
                    disp('Experiment aborted by user!');
                    disp(' ');                    
                    Screen('CloseAll');
                    save(fullfile(output,saveName)); % save everything
                    return;
            end
        end
    end
    
    % determine visibility next trial after every n trials
    if stairCount==nEvalTrials
        
        % trials to evaluate
        idx = (stairStep-1)*nEvalTrials+1:stairStep*nEvalTrials;
        
        V(stairStep) = vis_scale(visibility); % track
        acc(stairStep) = sum(trials(idx)==B(idx,trialResponse))/nEvalTrials;
        
        % update visibility based on accuracy
        if (acc(stairStep) > upper_bound) || (acc(stairStep) < lower_bound)
            visibility = visibility - round((acc(stairStep)-threshold)*120);
            if visibility > length(vis_scale); visibility = length(vis_scale); end % cap
        end
        
        % update counters
        stairStep = stairStep+1;
        stairCount = 0; % reset
    end    
    
    % Inter trial interval
    Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
        4, [255,255,255], [rect(3)/2, rect(4)/2], 1);
    Screen('Flip', w);
    WaitSecs(ITIs(iTrial));
    
    % Close all textures to free memory
    tmp = unique(cell2mat(target));
    for i_tex = 1:length(tmp)
        Screen('Close', tmp(i_tex));
    end
end
save(fullfile(output,saveName)); % save everything

Screen('TextSize',w, 28);
DrawFormattedText(w, 'This is the end of this calibration!', 'center', 'center', [255 255 255]);
vbl = Screen('Flip', w);
WaitSecs(2);
Screen('CloseAll')
sca; 