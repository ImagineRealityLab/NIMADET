function [PD,PA,V] = practiceDetection(subID,orientations)

%% practice detection (PD) script
% Opens window -> defines trials -> makes gabors -> adds noise ->
% runs through trials.

% =========================================================================
% Setup
% =========================================================================

% output
output = fullfile(cd,'results',subID);
if ~exist(output,'dir'); mkdir(output); end
saveName = sprintf('PD_%s.mat',subID);

[w, rect] = setWindow(0);
HideCursor;

% Trial number and P/A
nOri    = length(orientations);
nTrials = 8; % per orientation
PA      = [ones(nTrials/2,1); zeros(nTrials/2,1)];
PA      = repmat(PA,1,2);
PA(:,1) = PA(randperm(nTrials),1); PA(:,2) = PA(randperm(nTrials),2);
PD      = nan(nOri,nTrials,2);

% Visibility settings
threshold = 0.8;                        % acc to continue
vis_scale = [0 logspace(log10(0.005),log10(0.2),299)]; % steps in log space
visibility = 0.12;                      % Start visibility
[~,visibility] = min(abs(vis_scale-visibility));  % scale idx

% responses
trialResponse = 1;
trialRT       = 2;
yesKey        = 'h';
noKey         = 'j';

% timings    
fixTime       = 0.2;
mITI          = 2; % mean ITI - randomly sample from norm
sITI          = 1; % SD for sampling
ITIs          = normrnd(mITI,sITI,nOri,nTrials);

% =========================================================================
% Stimuli
% =========================================================================

% Makes the gabors to show for instruction
gaborPatch   = cell(nOri,1);
gaborTexture = cell(nOri,1);
for iOri = 1:nOri
    % stimulus
    gaborPatch{iOri} = make_stimulus(orientations(iOri),0.3);
    % texture
    gaborTexture{iOri} = Screen('MakeTexture',w,gaborPatch{iOri});
end

% Noise stimulus info
displayDuration = 2;                     % Duration of the stimulus in seconds
hz = Screen('NominalFrameRate', w);
ifi = 1/hz;                              % Refresh rate
nStepSize = 2;                           % 2 frames per step
nSteps = (displayDuration/ifi)/nStepSize;

cues = {'A','B'}; 

%% Instructon screen
[~, yCenter] = RectCenter(rect);
[x_pix, ~] = Screen('WindowSize', w);

% explain finger placement
text = ['During this experiment, you will have to use all your fingers'...
    ' except your thumbs to respond. \n \n Please place your left hand over the '...
    '[a],[s],[d],[f] keys, \n with your left pinky (LP) on the [a] and ' ...
    'left index (LI) on the [f]. \n \n Place your right hand over the [h],[j],[k],[l] '...
    ' \n with your right index (RI) on the [h] and your right pinky (RP) on the [l]. \n'...
    '\n [Press any key to continue] \n '];

Screen('TextSize',w, 28);
DrawFormattedText(w, text, 'center', yCenter*0.65, [255 255 255]);
vbl=Screen('Flip', w);
WaitSecs(1);
KbWait;

% show instructions
text = ['In this experiment you will be shown gratings in noise (see below) \n ',...
    'Your task is to indicate whether a grating is presented. \n ',...
    'Your right index finger (RI) indicates "Yes" and your right middle finger (RM) indicates "No" \n ',...
    'Please keep your eyes fixated on the fixation cross in the middle of the screen. \n ',...
    'A grating will be present in 50% of the trials. \n ',...
    'We will first practice this a few times. \n \n ',...
    '[Press any key to continue] \n '];

Screen('TextSize',w, 28);
DrawFormattedText(w, text, 'center', yCenter*0.65, [255 255 255]);

% show gratings
[xCenter, yCenter] = RectCenter(rect);
Xpos     = [x_pix*(1/3) x_pix*(2/3)];
baseRect = [0 0 size(gaborPatch{1},1) size(gaborPatch{1},1)];
allRects = nan(4, 3);
for i = 1:2
    allRects(:, i) = CenterRectOnPointd(baseRect, Xpos(i), yCenter*1.4);
end

Screen('DrawTextures', w, gaborTexture{1}, [], allRects(:,1), [],[], 0.5);
DrawFormattedText(w, 'Grating A', xCenter*(1.75/3), yCenter*1.125, [255 255 255]);

Screen('DrawTextures', w, gaborTexture{2}, [], allRects(:,2), [],[], 0.5);
DrawFormattedText(w, 'Grating B', xCenter*(3.75/3), yCenter*1.125, [255 255 255]);

vbl=Screen('Flip', w);
WaitSecs(1)
KbWait;

%% Trials start
for iOri = 1:nOri   
    notYet = 1;
    while notYet % until enough correct
    
        WaitSecs(fixTime); % otherwise it doesn't work
        
    % instruction screen with gratings
    text = sprintf('During this block, you will be detecting Grating %s (see below). \n After each trial, please indicate whether this grating was presented or not. \n \n [Press any key to start] \n ',cues{iOri});
    
    Screen('TextSize',w, 28);
    DrawFormattedText(w, text, 'center', yCenter*0.75, [255 255 255]);    
    
    Screen('DrawTextures', w, gaborTexture{1}, [], allRects(:,1), [],[], 0.5);
    DrawFormattedText(w, 'Grating A', xCenter*(1.75/3), yCenter*1.125, [255 255 255]);
    
    Screen('DrawTextures', w, gaborTexture{2}, [], allRects(:,2), [],[], 0.5);
    DrawFormattedText(w, 'Grating B', xCenter*(3.75/3), yCenter*1.125, [255 255 255]);
    
    vbl=Screen('Flip', w);
    KbWait;
    
    for iTrial = 1:nTrials
        
        % Fixation
        Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
            4, [0,0,0], [rect(3)/2, rect(4)/2], 1);
        vbl=Screen('Flip', w);
        WaitSecs(0.2);
        
        if PA(iTrial,iOri) == 1 % present trial
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
                make_stimulus(orientations(iOri),schedule(i_frame)));
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
                    PD(iOri,iTrial,trialResponse) = strcmp(key,yesKey); % 1 yes 0 no
                    PD(iOri,iTrial,trialRT) = keyTime-vbl;
                    
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
    end
    
    % check if the accuracy is high enough, otherwise go again
    acc = mean(squeeze(PD(iOri,:,trialResponse))'==PA(:,iOri));
    WaitSecs(0.2);
    if acc < threshold % too low
        visibility = visibility - round((acc-threshold)*100); % add visibility
        if visibility > length(vis_scale); visibility = length(vis_scale); end % cap at end
        text = sprintf('You were only correct on %d out of %d trials. \n So, let"s try this again! \n [Press any key to continue] \n ',sum(PD(iOri,:,trialResponse)'==PA(:,iOri)),nTrials);
        
        Screen('TextSize',w, 28);
        DrawFormattedText(w, text, 'center', 'center', 255);
        vbl = Screen('Flip', w);
        KbWait;
    else % good 
        notYet = 0; % we can continue
        text = sprintf('Well done! You were correct on %d out of %d trials. \n [Press any key to continue]',sum(PD(iOri,:,trialResponse)'==PA(:,iOri)),nTrials);
        %text = 'poop'; 
        Screen('TextSize',w, 28);
        DrawFormattedText(w, text, 'center', 'center', 255);
        vbl = Screen('Flip', w);
        KbWait;
    end   
    end
end

V = vis_scale(visibility); % use this for staircase and next practice 

save(fullfile(output,saveName)); % save everything

Screen('TextSize',w, 28);
DrawFormattedText(w, 'This is the end of the first practice task!', 'center', 'center', [255 255 255]);
vbl = Screen('Flip', w);
WaitSecs(2);
Screen('CloseAll')
sca;