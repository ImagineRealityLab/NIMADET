%% practice the main task
% contains all components and set at a relatively high visibility level,
% determined by the detection practice to lead to > 0.7 accuracy

function [R,C] = practiceMainTask(subID,orientations,vivDet,V)

% =========================================================================
% Setup
% =========================================================================
[w, rect] = setWindow(0);
HideCursor;

% output
output = fullfile(cd,'results',subID);
if ~exist(output,'dir'); mkdir(output); end
saveName = sprintf('PMT_%s.mat',subID);

% Trial numbers and order
nOri    = length(orientations);
nMB     = 2; % mini-blocks
nRep    = 1;
[blocks,miniblocks] = blockStructure(nOri,nMB,nRep);
nBlocks = numel(blocks);
nTrials = 8; % per mini-block
trials  = trialStructure(miniblocks,nTrials);

R       = nan(nBlocks*nMB,nTrials,4);
C       = nan(nBlocks,1); % ima check

% Visibility settings
vis_scale = [0 logspace(log10(0.005),log10(0.2),299)]; % steps in log space
[~,visibility] = min(abs(vis_scale-V));  % scale idx

% responses
vivRating     = 1;
vivRT         = 2;
detResponse   = 3;
detRT         = 4;
responseMappings = repmat(1:2,1,nBlocks/2);
responseMappings = responseMappings(randperm(nBlocks));

yesKey        = {'h','f'};
noKey         = {'j','d'};
vividnessKeys = {'a','s','d','f';'l','k','j','h'};
checkKeys     = {'k','l';'a','s'};

% timing
mITI    = 2; % mean ITI - randomly sample from norm
sITI    = 1; % SD for sampling
ITIs    = normrnd(mITI,sITI,nBlocks*nMB,nTrials);
fixTime = 0.2;

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

%% Instructon screen
[xCenter, yCenter] = RectCenter(rect);
[x_pix, ~] = Screen('WindowSize', w);

% show instructions 1.
text = ['Now for the main task, we will combine these two things together. \n ',...
    'During each block, you will imagine a grating of a specific orientation \n ',...
    'and indicate your vividness, like you did just now. \n \n',...
    'At the same time, you will be detecting a grating of a certain orientation, \n ',...
    'either the same as the one you are imagining, or the other one.\n \n ',...
    'You will imagine one orientation for an entire block, but the orientation \n ',...
    'you have to detect will switch halfway. We will tell you when this happens. \n \n '];

if vivDet == 1 % first vividness
    text = [text 'After each trial, you will first indicate how vivid your imagery was \n and then whether the to-be-detected grating was presented. \n \n'];
elseif vivDet == 2 % first detection
    text = [text 'After each trial, you will first indicate whether the to-be-detected grating was presented \n and then how vivid your imagery was. \n \n'];
end

text = [text '\n \n ',...
    '[Press any key to continue] \n '];

Screen('TextSize',w, 28);
DrawFormattedText(w, text, 'center', yCenter*0.75, [255 255 255]);

vbl=Screen('Flip', w);
WaitSecs(1);
KbWait;

% show instructions 2.

text = ['One other important change is that each block, which hand you will use \n'...
    'to indicate imagery vividness and which hand to indicate presence, switches. \n'...
    'You will be informed at the start of each block which hand to use for what. \n'...
    'Your index finger will always indicate yes and high vividness. \n \n '];

text = [text 'Now let"s practice all of this this a few times. \n \n ',...
    '[Press any key to continue] \n '];

Screen('TextSize',w, 28);
DrawFormattedText(w, text, 'center', yCenter*0.75, [255 255 255]);

vbl=Screen('Flip', w);
WaitSecs(1);
KbWait;

% grating info
% show gratings
[xCenter, yCenter] = RectCenter(rect);
Xpos     = [x_pix*(1/3) x_pix*(2/3)];
baseRect = [0 0 size(gaborPatch{1},1) size(gaborPatch{1},1)];
allRects = nan(4, 3);
for i = 1:2
    allRects(:, i) = CenterRectOnPointd(baseRect, Xpos(i), yCenter*1.4);
end
centerRect = CenterRectOnPointd(baseRect, xCenter, yCenter*1.4);

%% Trials start

miniblock_counter = 1;

for iBlock = 1:length(blocks)
    
    WaitSecs(fixTime);
    
    % instruction screen with gratings
    text = sprintf('This is block %d out of %d. \n \n',iBlock,nBlocks);
    text = [text sprintf('During this block, you will be imagining grating %s (see below). \n',cues{blocks(iBlock)})];
    text = [text 'Please imagine this grating as vividly as possible during each trial, \n as if it was actually presented. \n \n'];
    
    if responseMappings(iBlock) == 1
        text = [text 'You will use your left hand to indicate vividness \n and your right to indicate whether a stimulus was presented.'];
        RM = 1;
    elseif responseMappings(iBlock) == 2
        text = [text 'You will use your right hand to indicate vividness \n and your left to indicate whether a stimulus was presented.'];
        RM = 2;
    end
    text = [text '\n \n [Press any key to start] \n '];
    
    Screen('TextSize',w, 28);
    DrawFormattedText(w, text, 'center', yCenter*0.6, [255 255 255]);
    
    Screen('DrawTextures', w, gaborTexture{blocks(iBlock)}, [], centerRect, [],[], 0.5);
    DrawFormattedText(w, sprintf('Grating %s',cues{blocks(iBlock)}), xCenter*0.95, yCenter*1.2, [255 255 255]);
    
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
                    
                    vbl=Screen('Flip', w);
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
                vbl = Screen('Flip', w);
                
                keyPressed = 0; % clear previous response
                while ~keyPressed
                    
                    [~, keyTime, keyCode] = KbCheck(-3);
                    key = KbName(keyCode);
                    
                    if ~iscell(key) % only start a keypress if there is only one key being pressed
                        if any(strcmp(key, vividnessKeys(RM,:)))
                            
                            % fill in B
                            R(miniblock_counter,iTrial,vivRating) = find(strcmp(key,vividnessKeys(RM,:))); % 1 to 4
                            R(miniblock_counter,iTrial,vivRT) = keyTime-vbl; % RT
                            
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
                vbl = Screen('Flip', w);
                
                % Log response
                keyPressed = 0; % clear previous response
                while ~keyPressed
                    
                    [~, keyTime, keyCode] = KbCheck(-3);
                    key = KbName(keyCode);
                    
                    if ~iscell(key) % only start a keypress if there is only one key being pressed
                        if any(strcmp(key, {yesKey{RM},noKey{RM}}))
                            
                            % fill in B
                            R(miniblock_counter,iTrial,detResponse) = strcmp(key,yesKey{RM}); % 1 yes 0 no
                            R(miniblock_counter,iTrial,detRT) = keyTime-vbl;
                            
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
                vbl = Screen('Flip', w);
                
                % Log response
                keyPressed = 0; % clear previous response
                while ~keyPressed
                    
                    [~, keyTime, keyCode] = KbCheck(-3);
                    key = KbName(keyCode);
                    
                    if ~iscell(key) % only start a keypress if there is only one key being pressed
                        if any(strcmp(key, {yesKey{RM},noKey{RM}}))
                            
                            % fill in B
                            R(miniblock_counter,iTrial,detResponse) = strcmp(key,yesKey{RM}); % 1 yes 0 no
                            R(miniblock_counter,iTrial,detRT) = keyTime-vbl;
                            
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
                vbl = Screen('Flip', w);
                
                keyPressed = 0; % clear previous response
                while ~keyPressed
                    
                    [~, keyTime, keyCode] = KbCheck(-3);
                    key = KbName(keyCode);
                    
                    if ~iscell(key) % only start a keypress if there is only one key being pressed
                        if any(strcmp(key, vividnessKeys(RM,:)))
                            
                            % fill in B
                            R(miniblock_counter,iTrial,vivRating) = find(strcmp(key,vividnessKeys(RM,:))); % 1 to 4
                            R(miniblock_counter,iTrial,vivRT) = keyTime-vbl; % RT
                            
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
        
        % update miniblock counter
        miniblock_counter = miniblock_counter+1;
    end
    
    % Imagery check
    text = 'CHECK! Which grating did you imagine this block? \n';
    if RM == 1 % Right hand
        text = [text 'Grating A [RR] or Grating B [RP]'];
    else
        text = [text 'Grating A [LP] or Grating B [LR]'];
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
        
        [~, keyTime, keyCode] = KbCheck(-3);
        key = KbName(keyCode);
        
        if ~iscell(key) % only start a keypress if there is only one key being pressed
            if any(strcmp(key, checkKeys(RM,:)))
                
                % fill in response
                checkResponse = find(strcmp(key,checkKeys(RM,:))); % 1 to 2
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
        text = 'Correct! \n \n [Press any key to continue]';
    elseif C(iBlock) == 0
        text = 'That is incorrect, please read the block instructions carefully! \n \n [Press any key to continue]';
    end
    Screen('TextSize',w, 28);
    DrawFormattedText(w, text, 'center', 'center', 255);
    vbl = Screen('Flip', w);
    KbWait;
end

save(fullfile(output,saveName)); % save everything


Screen('TextSize',w, 28);
DrawFormattedText(w, 'This is the end of the final practice task!', 'center', 'center', [255 255 255]);
vbl = Screen('Flip', w);
WaitSecs(2);
Screen('CloseAll')
sca;
