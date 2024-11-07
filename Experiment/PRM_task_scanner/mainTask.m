%% main task script

function [R,C,T] = mainTask(subID,orientations,visibility,vivDet,environment,run,RM)

% =========================================================================
% Setup
% =========================================================================
[w, rect] = setWindow(0);
HideCursor;

% output
output = fullfile(cd,'results',subID);
if ~exist(output,'dir'); mkdir(output); end
saveName = sprintf('MT_%s_run%d.mat',subID,run);

% Trials and blocks
nOri    = length(orientations);
nMB     = 2; % mini-blocks
nRep    = 2;
[blocks,miniblocks] = blockStructure(nOri,nMB,nRep);
nBlocks = numel(blocks);
nTrials = 12; % per mini-block
trials  = trialStructure(miniblocks,nTrials);

% Visibility settings
vis_scale = [0 logspace(log10(0.005),log10(0.2),299)]; % steps in log space
[~,visibility(1)] = min(abs(vis_scale-visibility(1)));  % scale idx
[~,visibility(2)] = min(abs(vis_scale-visibility(2)));

% responses
R       = nan(nBlocks*nMB,nTrials,4);
C       = nan(nBlocks,1); % ima check

% timings
mITI    = 2; % mean ITI - randomly sample from norm
sITI    = 1; % SD for sampling
ITIs    = normrnd(mITI,sITI,nBlocks*nMB,nTrials);

cStim   = 1; % timing of stim onset
cDet    = 2; % detection screen
cViv    = 3; % vividness screen
T.presTimes = nan(nBlocks,nTrials,3);

fixTime = 0.2;
checkfeedbackTime = 2; % how long to show the imagery check feedback

% response keys
vivRating     = 1;
vivRT         = 2;
detResponse   = 3;
detRT         = 4;
if strcmp(environment,'mri')
    if RM == 1
        yesKey        = '1!';
        noKey         = '2@';
        vividnessKeys = {'9(','8*','7&','6^'}; % left hand
        checkKeys     = {'3#','4$'};
    elseif RM == 2
        yesKey        = '6^';
        noKey         = '7&';
        vividnessKeys = {'1!','2@','3#','4$'}; % right hand
        checkKeys     = {'9(','8*'};
    end
else
    if RM == 1
        yesKey        = 'h';
        noKey         = 'j';
        vividnessKeys = {'a','s','d','f'};
        checkKeys     = {'k','l'};
    elseif RM == 2
        yesKey        = 'f';
        noKey         = 'd';
        vividnessKeys = {'l','k','j','h'};
        checkKeys     = {'a','s'};
    end
end

% only check for certain keys
RestrictKeys = {yesKey,noKey,vividnessKeys{:},checkKeys{:},'5%','SPACE','ESCAPE'};
RestrictKeyCodes = nan(length(RestrictKeys),1);
for k = 1:length(RestrictKeys); RestrictKeyCodes(k) = KbName(RestrictKeys{k}); end
RestrictKeysForKbCheck(RestrictKeyCodes);

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

cues = {'A','B'};


%% Wait for scanner trigger
[xCenter, yCenter] = RectCenter(rect);
[x_pix, ~] = Screen('WindowSize', w);

triggerTimeSecs = GetSecs;
KbQueueRelease; KbQueueRelease; % clear queue prior to start
ListenChar(0);    

if strcmp(environment,'mri') % wait for scanner trigger
    excludeVolumes = 5; % slicesPerVolume = 48;
    scannerTrigger = KbName('5%');
    
    num_five = 0; % how many triggers did I get from the scanner?
    while num_five < excludeVolumes % * slicesPerVolume
        
        text = 'Please wait for the scanner to start...';
        DrawFormattedText(w, text, 'center', 'center', [255 255 255]);
        vbl = Screen('Flip', w);
        
        triggerTimeSecs = KbTriggerWait(scannerTrigger);
        num_five = num_five + 1;
        disp(['Detected volume ' num2str(num_five)])
    end
    
    % now don't listen to the scanner trigger anymore!
    DisableKeysForKbCheck(scannerTrigger)
    
end
ListenChar(1);  

% clear queue again
%KbQueueRelease; KbQueueRelease;

% log start time
T.starttime = triggerTimeSecs;
WaitSecs(fixTime);

%% Show instruction screen

text = sprintf('This is run number %d \n',run);
if RM == 1
    text = [text 'During this whole run, you will use your left hand to indicate vividness \n and your right to indicate whether a stimulus was presented. \n \n'];
elseif RM == 2
    text = [text 'During this whole run, you will use your right hand to indicate vividness \n and your left to indicate whether a stimulus was presented. \n \n'];
end
text = [text '[Press any key to continue]'];
Screen('TextSize',w, 28);
DrawFormattedText(w, text, 'center', 'center', [255 255 255]);
vbl = Screen('Flip', w);
KbWait;

% info to show gratings
Xpos     = [x_pix*(1/3) x_pix*(2/3)];
baseRect = [0 0 size(gaborPatch{1},1) size(gaborPatch{1},1)];

allRects = nan(4, 3);
for i = 1:2
    allRects(:, i) = CenterRectOnPointd(baseRect, Xpos(i), yCenter*1.4);
end


%% Trials start

miniblock_counter = 1;

for iBlock = 1:length(blocks)
    
    WaitSecs(fixTime);
    
    % instruction screen with gratings
    text = sprintf('This is block %d out of %d of this run. \n \n',iBlock,nBlocks);
    text = [text sprintf('During this block, you will be imagining grating %s (see below). \n',cues{blocks(iBlock)})];
    text = [text 'Please remember to imagine this grating as vividly as possible during each trial, \n as if it was actually presented, \n'];
    text = [text 'and to keep your eyes fixated on the cross in the center of the screen. \n \n'];
    text = [text '[Press any key to continue] \n '];
    
    Screen('TextSize',w, 28);
    DrawFormattedText(w, text, 'center', yCenter*0.75, [255 255 255]);
    
    Screen('DrawTextures', w, gaborTexture{1}, [], allRects(:,1), [],[], 0.5);
    DrawFormattedText(w, 'Grating A', xCenter*(1.8/3), yCenter*1.2, [255 255 255]);
    
    Screen('DrawTextures', w, gaborTexture{2}, [], allRects(:,2), [],[], 0.5);
    DrawFormattedText(w, 'Grating B', xCenter*(3.8/3), yCenter*1.2, [255 255 255]);
    
    vbl=Screen('Flip', w);
    WaitSecs(1);
    KbWait;
    
    % loop over miniblocks
    for iMiniBlock = 1:nMB
        
        text = sprintf('During the next few trials you will be detecting grating %s \n',cues{miniblocks(miniblock_counter)});
        text = [text sprintf('Remember to still imagine grating %s during each trial. \n \n',cues{blocks(iBlock)})];
        text = [text 'Remember that a grating will be present in 50% of the trials. \n You will need to concentrate hard as the gratings will be difficult to see \n'];
        
        text = [text '\n \n [Press any key to start] \n '];
        
        Screen('TextSize',w, 28);
        DrawFormattedText(w, text, 'center', yCenter*0.75, [255 255 255]);
        
        Screen('DrawTextures', w, gaborTexture{1}, [], allRects(:,1), [],[], 0.5);
        DrawFormattedText(w, 'Grating A', xCenter*(1.75/3), yCenter*1.125, [255 255 255]);
        
        Screen('DrawTextures', w, gaborTexture{2}, [], allRects(:,2), [],[], 0.5);
        DrawFormattedText(w, 'Grating B', xCenter*(3.75/3), yCenter*1.125, [255 255 255]);
        
        vbl=Screen('Flip', w);
        WaitSecs(1);
        KbWait;
        
        % loop over trials
        for iTrial = 1:nTrials
            
            % Fixation
            Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
                4, [0,0,0], [rect(3)/2, rect(4)/2], 1);
            vbl=Screen('Flip', w);
            WaitSecs(fixTime);
            
            if trials(miniblock_counter,iTrial) == 1 % present trial
                % schedule of visibility gradient (i.e how visible at each frame)
                % Increases till most visible at the end
                schedule = vis_scale(round(linspace(1,visibility(1),nSteps)));
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
                    make_stimulus(orientations(miniblocks(miniblock_counter,1)),schedule(i_frame)));
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
                    
                    if i_frame == 1 % log start-time
                        T.presTimes(miniblock_counter,iTrial,cStim) = Screen('Flip', w);
                    else
                        vbl=Screen('Flip', w);
                    end
                end
            end
            
            if vivDet == 1
                
                % Vividness rating first
                text = 'How vivid was your imagery? \n';
                if RM == 2 % right hand
                    text = [text '\n 4[RI] - 1[RL]'];
                else
                    text = [text '\n 1[LL] - 4[LI]'];
                end
                Screen('TextSize',w, 28);
                DrawFormattedText(w, text, 'center', 'center', 255);
                T.presTimes(miniblock_counter,iTrial,cViv) = Screen('Flip', w);
                
                keyPressed = 0; % clear previous response
                while ~keyPressed
                    
                    [~, keyTime, keyCode] = KbCheck(-3);
                    key = KbName(keyCode);
                    
                    if ~iscell(key) % only start a keypress if there is only one key being pressed
                        if any(strcmp(key, vividnessKeys))
                            
                            % fill in B
                            R(miniblock_counter,iTrial,vivRating) = find(strcmp(key,vividnessKeys)); % 1 to 4
                            R(miniblock_counter,iTrial,vivRT) = keyTime-T.presTimes(miniblock_counter,iTrial,cViv); % RT
                            
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
                
                % Detection
                text = 'Was there a grating on the screen? \n';
                if RM == 1 % right hand
                    text = [text 'Yes [RI] or no [RM]'];
                else
                    text = [text 'No [LM] or yes [LI]'];
                end
                Screen('TextSize',w, 28);
                DrawFormattedText(w, text, 'center', 'center', 255);
                T.presTimes(miniblock_counter,iTrial,cDet) = Screen('Flip', w);
                
                % Log response
                keyPressed = 0; % clear previous response
                while ~keyPressed
                    
                    [~, keyTime, keyCode] = KbCheck(-3);
                    key = KbName(keyCode);
                    
                    if ~iscell(key) % only start a keypress if there is only one key being pressed
                        if any(strcmp(key, {yesKey,noKey}))
                            
                            % fill in B
                            R(miniblock_counter,iTrial,detResponse) = strcmp(key,yesKey); % 1 yes 0 no
                            R(miniblock_counter,iTrial,detRT) = keyTime-T.presTimes(miniblock_counter,iTrial,cDet);
                            
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
                
            elseif vivDet == 2
                
                % Detection first
                text = 'Was there a grating on the screen? \n';
                if RM == 1 % right hand
                    text = [text 'Yes [RI] or no [RM]'];
                else
                    text = [text 'No [LM] or yes [LI]'];
                end
                Screen('TextSize',w, 28);
                DrawFormattedText(w, text, 'center', 'center', 255);
                T.presTimes(miniblock_counter,iTrial,cDet) = Screen('Flip', w);
                
                % Log response
                keyPressed = 0; % clear previous response
                while ~keyPressed
                    
                    [~, keyTime, keyCode] = KbCheck(-3);
                    key = KbName(keyCode);
                    
                    if ~iscell(key) % only start a keypress if there is only one key being pressed
                        if any(strcmp(key, {yesKey,noKey}))
                            
                            % fill in B
                            R(miniblock_counter,iTrial,detResponse) = strcmp(key,yesKey); % 1 yes 0 no
                            R(miniblock_counter,iTrial,detRT) = keyTime-T.presTimes(miniblock_counter,iTrial,cDet);
                            
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
                
                % Vividness rating
                text = 'How vivid was your imagery? \n';
                if RM == 2 % right hand
                    text = [text '\n 4[RI] - 1[RL]'];
                else
                    text = [text '\n 1[LL] - 4[LI]'];
                end
                Screen('TextSize',w, 28);
                DrawFormattedText(w, text, 'center', 'center', 255);
                T.presTimes(miniblock_counter,iTrial,cViv) = Screen('Flip', w);
                
                keyPressed = 0; % clear previous response
                while ~keyPressed
                    
                    [~, keyTime, keyCode] = KbCheck(-3);
                    key = KbName(keyCode);
                    
                    if ~iscell(key) % only start a keypress if there is only one key being pressed
                        if any(strcmp(key, vividnessKeys))
                            
                            % fill in B
                            R(miniblock_counter,iTrial,vivRating) = find(strcmp(key,vividnessKeys)); % 1 to 4
                            R(miniblock_counter,iTrial,vivRT) = keyTime-T.presTimes(miniblock_counter,iTrial,cViv); % RT
                            
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
            end
            
            % Inter trial interval
            Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
                4, [255,255,255], [rect(3)/2, rect(4)/2], 1);
            Screen('Flip', w);
            WaitSecs(ITIs(miniblock_counter,iTrial));
            
            % Close all textures to free memory
            tmp = unique(cell2mat(target));
            for i_tex = 1:length(tmp)
                Screen('Close', tmp(i_tex));
            end
            
        end
        
        % update minminiblock_counter counter
        miniblock_counter = miniblock_counter+1;
        
    end
    
    % Imagery check
    text = 'CHECK! Which grating did you imagine this block? \n';
    if RM == 1 % Right hand
        text = [text 'Grating A [RR] or Grating B [RL]'];
    else
        text = [text 'Grating A [LL] or Grating B [LR]'];
    end
    Screen('TextSize',w, 28);
    DrawFormattedText(w, text, 'center', yCenter, [255 255 255]);
    
    Screen('DrawTextures', w, gaborTexture{1}, [], allRects(:,1), [],[], 0.5);
    DrawFormattedText(w, 'Grating A', xCenter*(1.8/3), yCenter*1.2, [255 255 255]);
    
    Screen('DrawTextures', w, gaborTexture{2}, [], allRects(:,2), [],[], 0.5);
    DrawFormattedText(w, 'Grating B', xCenter*(3.8/3), yCenter*1.2, [255 255 255]);
    
    vbl=Screen('Flip', w);
    
    % log response
    keyPressed = 0; % clear previous response
    while ~keyPressed
        
        [~, ~, keyCode] = KbCheck(-3);
        key = KbName(keyCode);
        
        if ~iscell(key) % only start a keypress if there is only one key being pressed
            if any(strcmp(key, checkKeys))
                
                % fill in response
                checkResponse = find(strcmp(key,checkKeys)); % 1 to 2
                if checkResponse == blocks(iBlock) % correct
                    C(iBlock) = 1;
                else
                    C(iBlock) = 0;
                end
                
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
    WaitSecs(0.2);
    if C(iBlock) == 1
        text = 'Correct!';
    elseif C(iBlock) == 0
        text = 'That is incorrect, please read the block instructions carefully!';
    end
    Screen('TextSize',w, 28);
    DrawFormattedText(w, text, 'center', 'center', 255);
    vbl = Screen('Flip', w);
    WaitSecs(checkfeedbackTime);
    
    % Break
    text = sprintf('This is the end of block %d of %d for this run. \n You can now have a short break. \n \n \n [Press any key when you are ready to continue] \n ',iBlock,nBlocks);
    
    Screen('TextSize',w, 28);
    DrawFormattedText(w, text, 'center', 'center', [255 255 255]);
    
    vbl=Screen('Flip', w);
    KbWait;
    
end

save(fullfile(output,saveName)); % save everything


Screen('TextSize',w, 28);
DrawFormattedText(w, 'This is the end of this run. Please relax while we restart the scanner...', 'center', 'center', [255 255 255]);
vbl = Screen('Flip', w);
WaitSecs(2);
Screen('CloseAll');
sca;
