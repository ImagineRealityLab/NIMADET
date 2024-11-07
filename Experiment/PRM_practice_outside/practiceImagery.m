%% practice the imagery task
% show which cue goes with which grating
% only catch trials - train people to do the imagery task well

function [R_V,R_P] = practiceImagery(subID,orientations)

% =========================================================================
% Setup
% =========================================================================
[w, rect] = setWindow(0);
HideCursor;

% output
output = fullfile(cd,'results',subID);
if ~exist(output,'dir'); mkdir(output); end
saveName = sprintf('PI_%s.mat',subID);

% Trial number and percentage of catch trials
nOri        = length(orientations);
nTrialsC    = 5; % how many correct in a row?
nRemind     = 6; % show gratings every N trials
nTrials     = 100;
offsets     = nan(nOri,nTrials);
tmp         = normrnd(0,10,nTrials*100,1); tmp = tmp(abs(tmp) > 5.5); % otherwise it is too small
offsets(1,:) = tmp(1:nTrials); offsets(2,:) = tmp(end-nTrials+1:end);

% cues
cues = {'A','B'};

% keys
ccKey = 'h';
cwKey = 'j';
vividnessKeys = {'a','s','d','f'};

% responses
R_V = zeros(nOri,nTrials,2);
R_P = zeros(nOri,nTrials,3);

% timing
probeTime = 0.5;
fixTime   = 0.2;
mITI    = 2; % mean ITI - randomly sample from norm
sITI    = 1; % SD for sampling
ITIs    = normrnd(mITI,sITI,nOri,nTrials);

% =========================================================================
% Stimuli
% =========================================================================

% Makes the gabors to show for instruction
gaborPatch   = cell(nOri,1);
gaborTexture = cell(nOri,1);
for iOri = 1:nOri
    % stimulus
    gaborPatch{iOri} = make_stimulus(orientations(iOri),1); % full visibility
    % texture
    gaborTexture{iOri} = Screen('MakeTexture',w,gaborPatch{iOri});
end

% Noise stimulus info
displayDuration = 2;                     % Duration of the stimulus in seconds
hz = Screen('NominalFrameRate', w);
ifi = 1/hz;                              % Refresh rate
nStepSize = 2;                           % 2 frames per step
nSteps = (displayDuration/ifi)/nStepSize;

%% Instructon screen
[xCenter, yCenter] = RectCenter(rect);
[x_pix, ~] = Screen('WindowSize', w);

% show instructions (screen 1)
text = ['In this part of the experiment you will imagine the gratings \n ',...
    'while again looking at dynamicly changing noise (so keep your eyes open!). \n ',...
    'On each trial, you should start imagining when the fixation cross turns black. \n \n',...
    'Your task is to first indicate how vivid your mental image was\n ',...
    'on a scale from 1-4, using your left hand (the [a],[s],[d],[f] keys). \n  \n ',...
    '[Press any key to continue] \n '];

Screen('TextSize',w, 28);
DrawFormattedText(w, text, 'center', yCenter*0.6, [255 255 255]);

% show gratings
Xpos     = [x_pix*(1/3) x_pix*(2/3)];
baseRect = [0 0 size(gaborPatch{1},1) size(gaborPatch{1},1)];

allRects = nan(4, 3);
for i = 1:2
    allRects(:, i) = CenterRectOnPointd(baseRect, Xpos(i), yCenter*1.4);
end

Screen('DrawTextures', w, gaborTexture{1}, [], allRects(:,1), [],[], 0.5);
DrawFormattedText(w, 'Grating A', xCenter*(1.8/3), yCenter*1.2, [255 255 255]);

Screen('DrawTextures', w, gaborTexture{2}, [], allRects(:,2), [],[], 0.5);
DrawFormattedText(w, 'Grating B', xCenter*(3.8/3), yCenter*1.2, [255 255 255]);

vbl=Screen('Flip', w);
WaitSecs(1);
KbWait;

% show instructions (screen 2)
text = ['After that, another grating will be presented which is tilted with \n'...
    'respect to the grating you just imagined. \n \n ',...
    'You have to indicate whether the tilt is counter-clockwise (leftward) [RI] \n'...
    'or clockwise (rightward) [RM] with respect to the imagined grating. \n \n',...
    'You will practice imagining each orientation until you have 5 correct responses in a row. \n \n ',...
    '[Press any key to continue] \n '];

Screen('TextSize',w, 28);
DrawFormattedText(w, text, 'center', yCenter*0.6, [255 255 255]);

% show gratings
Xpos     = [x_pix*(1/3) x_pix*(2/3)];
baseRect = [0 0 size(gaborPatch{1},1) size(gaborPatch{1},1)];

allRects = nan(4, 3);
for i = 1:2
    allRects(:, i) = CenterRectOnPointd(baseRect, Xpos(i), yCenter*1.4);
end

Screen('DrawTextures', w, gaborTexture{1}, [], allRects(:,1), [],[], 0.5);
DrawFormattedText(w, 'Grating A', xCenter*(1.8/3), yCenter*1.2, [255 255 255]);

Screen('DrawTextures', w, gaborTexture{2}, [], allRects(:,2), [],[], 0.5);
DrawFormattedText(w, 'Grating B', xCenter*(3.8/3), yCenter*1.2, [255 255 255]);

vbl=Screen('Flip', w);
WaitSecs(1);
KbWait;
%% Trials start
for iOri = 1:nOri
    
    WaitSecs(fixTime);
    
    % instruction screen with gratings
    text = sprintf('During the next trials, please imagine Grating %s (see below) \n as vividly as possible, as if it was really appearing in the noise. \n \n Keep your eyes fixated on the fixation cross as much as possible. \n \n [Press any key to start] \n ',cues{iOri});
    
    Screen('TextSize',w, 28);
    DrawFormattedText(w, text, 'center', yCenter*0.75, [255 255 255]);    
    
    Screen('DrawTextures', w, gaborTexture{1}, [], allRects(:,1), [],[], 0.5);
    DrawFormattedText(w, 'Grating A', xCenter*(1.8/3), yCenter*1.2, [255 255 255]);
    
    Screen('DrawTextures', w, gaborTexture{2}, [], allRects(:,2), [],[], 0.5);
    DrawFormattedText(w, 'Grating B', xCenter*(3.8/3), yCenter*1.2, [255 255 255]);
    
    vbl=Screen('Flip', w);
    WaitSecs(1);
    KbWait;
    
    count = 0; iTrial = 0;
    while count < nTrialsC && iTrial < nTrials
        iTrial = iTrial + 1;
        
        % Fixation
        Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
            4, [0,0,0], [rect(3)/2, rect(4)/2], 1);
        vbl=Screen('Flip', w);
        WaitSecs(fixTime);
        
        % Make the textures for dynamic noise
        schedule = zeros(1,nSteps);
        target = {}; 
        for i_frame = 1:nSteps
            idx = ((i_frame-1)*nStepSize)+1:(i_frame*nStepSize);
            tmp = Screen('MakeTexture',w, ...
                make_stimulus(orientations(iOri),schedule(i_frame)));
            for i = 1:length(idx); target{idx(i)} = tmp; end
        end
        
        % Present dynamic noise
        tini = GetSecs;
        for i_frame = 1:length(target)
            
            while GetSecs-tini < ifi*i_frame
                
                Screen('DrawTextures',w,target{i_frame});
                
                Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
                    4, [0,0,0], [rect(3)/2, rect(4)/2], 1);
                
                vbl=Screen('Flip', w);
            end
        end
        
        % Vividness rating
        text = 'How vivid was your imagery? \n 1[LL] - 4[LI]';
        Screen('TextSize',w, 28);
        DrawFormattedText(w, text, 'center', 'center', 255);
        vbl = Screen('Flip', w);
        
        keyPressed = 0; % clear previous response
        while ~keyPressed
            
            [~, keyTime, keyCode] = KbCheck(-3);
            key = KbName(keyCode);
            
            if ~iscell(key) % only start a keypress if there is only one key being pressed
                if any(strcmp(key, vividnessKeys))
                    
                    % fill in B
                    R_V(iOri,iTrial,1) = find(strcmp(key,vividnessKeys)); % 1 to 4
                    R_V(iOri,iTrial,2) = keyTime-vbl; % RT
                    
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
        
        % Present probe
        Screen('DrawTextures', w, gaborTexture{iOri}, [], [], offsets(iOri,iTrial),[], 0.5);
        
        Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
            4, [0,0,0], [rect(3)/2, rect(4)/2], 1);
        
        vbl=Screen('Flip', w);
        WaitSecs(probeTime);
        
        % Probe discrimination
        text = 'Counter-clockwise [RI] or clockwise [RM]?';
        Screen('TextSize',w, 28);
        DrawFormattedText(w, text, 'center', 'center', 255);
        vbl = Screen('Flip', w);
        
        keyPressed = 0; % clear previous response
        while ~keyPressed
            
            [~, keyTime, keyCode] = KbCheck(-3);
            key = KbName(keyCode);
            
            if ~iscell(key) % only start a keypress if there is only one key being pressed
                if any(strcmp(key, {cwKey,ccKey}))
                    
                    % fill in B
                    R_P(iOri,iTrial,1) = find(strcmp(key,{cwKey,ccKey})); % 1 cc 2 cw
                    R_P(iOri,iTrial,2) = keyTime-vbl;
                    
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
        
        % Feedback
        if (offsets(iOri,iTrial) > 0) == (R_P(iOri,iTrial,1)==1)
            R_P(iOri,iTrial,3) = 1;
            count = count+1;
            text = sprintf('Correct! %d in a row!',count);
        elseif (offsets(iOri,iTrial) > 0) ~= (R_P(iOri,iTrial,1)==1)
            text = 'Incorrect';
            count = 0;
        end
        Screen('TextSize',w, 28);
        DrawFormattedText(w, text, 'center', 'center', 255);
        vbl = Screen('Flip', w);
        WaitSecs(0.5)
        
        % Inter trial interval 
        Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
            4, [255,255,255], [rect(3)/2, rect(4)/2], 1);
        Screen('Flip', w);
        WaitSecs(ITIs(iOri,iTrial));  
        
        % Close all textures to free memory
        tmp = unique(cell2mat(target));
        for i_tex = 1:length(tmp)
            Screen('Close', tmp(i_tex));
        end
        
        % show gratings as a reminder
        if mod(iTrial,nRemind) == 0
            
            text = sprintf('Just a reminder of what the gratings look like! \n You should imagine Grating %s during this block. \n \n [Press any key to continue] \n ',cues{iOri});
    
            Screen('TextSize',w, 28);
            DrawFormattedText(w, text, 'center', yCenter*0.75, [255 255 255]);
            
            Screen('DrawTextures', w, gaborTexture{1}, [], allRects(:,1), [],[], 0.5);
            DrawFormattedText(w, 'Grating A', xCenter*(1.8/3), yCenter*1.2, [255 255 255]);
            
            Screen('DrawTextures', w, gaborTexture{2}, [], allRects(:,2), [],[], 0.5);
            DrawFormattedText(w, 'Grating B', xCenter*(3.8/3), yCenter*1.2, [255 255 255]);
            
            vbl=Screen('Flip', w);
            KbWait;
        end
        
    end
end

save(fullfile(output,saveName)); % save everything

Screen('TextSize',w, 28);
DrawFormattedText(w, 'This is the end of the imagery practice!', 'center', 'center', [255 255 255]);
vbl = Screen('Flip', w);
WaitSecs(2);
Screen('CloseAll')
sca;