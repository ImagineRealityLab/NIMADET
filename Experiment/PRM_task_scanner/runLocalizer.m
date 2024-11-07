function runLocalizer(subID,orientations,TR,environment,run)

% =========================================================================
% Setup
% =========================================================================
[w, rect] = setWindow(0);
ShowCursor;

% output
output = fullfile(cd,'Results',subID);
saveName = sprintf('localizer_%s_run%d.mat',subID,run);

% trial numbers and timing
nOri    = length(orientations);
nBlocks = 7; % repetitions per orientation  - ~1.33 min per rep
seqOri  = nOri*2; % orientation1 - blank - orientation2 - blank - orientation3 - blank
presTime = 9*TR; % how many TR's to present one stimulus at a time 
flickerFreq = 4;
mBlockS = ceil(presTime*flickerFreq); % how long each mini block lasts in seconds
stimOrder = nan(nBlocks,seqOri);
T.presentationTime = zeros(nBlocks,seqOri,mBlockS);

% responses
if strcmp(environment,'mri')
    resKey = '1!'; % right index
else
    resKey = 'h';
end
R      = nan(nBlocks,4,floor(presTime*flickerFreq));

% =========================================================================
% Stimuli
% =========================================================================

% Predraw the gabors with different phases
gaborPatch   = cell(nOri,2); % one normal and one blue
gaborTexture = cell(nOri,2);
colours      = [0 3]; % not coloured and blue
for iOri = 1:nOri
    for c = 1:2
        % stimulus
        gaborPatch{iOri,c} = make_localizer_stimulus(orientations(iOri),colours(c)); % full visibility
        % texture
        gaborTexture{iOri,c} = Screen('MakeTexture',w,gaborPatch{iOri,c});
    end
end

%% Instructon screen + wait for scanner trigger
%instruction task
if run == 1
text = 'This the last part of the experiment. \n During this part, you will again be looking at oriented gratings. \n However, this time they are not embedded in noise. \n Instead of indicating whether a grating was present or not, \n now your task is to respond whenever the grating turns blue. \n When it turns blue, please press the right index button as fast as possible. \n Keep your eyes fixated at the cross during the entire task. \n This is run 1 out of 2 \n \n [Press any key to continue]';
else
    text = 'This is run 2/2 \n \n [Press any key to continue]';
end

Screen('TextSize',w, 28);
DrawFormattedText(w, text, 'center', 'center', [255 255 255]);
vbl=Screen('Flip', w);
KbWait;
WaitSecs(0.2);

% wait for trigger
excludeVolumes = 5; 
scannerTrigger = KbName('5%');

%initialize
num_five = 0; % how many triggers did I get from the scanner?
KbQueueRelease; KbQueueRelease; % clear cue prior to start
ListenChar(0); 
while num_five < excludeVolumes % * slicesPerVolume
    
    text = 'Please wait for the scanner to start... \n ';
    DrawFormattedText(w, text, 'center', 'center', [255 255 255]);
    vbl = Screen('Flip', w);    
    
    triggerTimeSecs = KbTriggerWait(scannerTrigger);
    num_five = num_five + 1;
    disp(['Detected volume ' num2str(num_five)]) 
end  
ListenChar(1); 

% log start time
T.starttime = triggerTimeSecs;
time = T.starttime;

%% Start localizer
% Remove any keypresses that occured before presentation of the
% stimuli.
FlushEvents('keyDown');
responseCount = 0;

changeTimesAbs = [];
changeTimes    = cell(nBlocks,1);

% generate random sequence of stimulus presentation
for iBlock = 1:nBlocks
    
    oriOrder = randperm(nOri);
    orientationOrder    = [oriOrder(1,1) NaN oriOrder(1,2) NaN oriOrder(1,3) NaN];
    stimOrder(iBlock,:) = orientationOrder;
    
    % How many times will the grating change colour in this
    % block? between 5-10 times.
    nChanges = 5 + round(5*rand);
    % When will these changes occur?
    changeTimes{iBlock} = randperm(floor(presTime*4-2));
    changeTimes{iBlock} = sort(changeTimes{iBlock}(1:nChanges));
    
    % Present the four min-blocks of 32s each: the two orientations interleaved with
    % fixation blocks.
    for i = 1:seqOri
        
        % Decide at which timepoint (in seconds) the fixation point will turn
        % black
        taskPres = changeTimes{iBlock}(changeTimes{iBlock} >= presTime*(i-1) & changeTimes{iBlock} <= presTime*i);
        % During which stimulus presentation (i.e. flicker)?
        taskPres = round((taskPres - presTime*(i-1)) * flickerFreq);
        
        for j = 1:mBlockS
                        
            if mod(j,2) && ~isnan(orientationOrder(i)) % when present grating
                
                Screen('DrawTextures',w,gaborTexture{stimOrder(iBlock,i),1});
                Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
                    4, [0,0,0], [rect(3)/2, rect(4)/2], 1);
                
            end
            if find(taskPres == j) % blue                
                if ~isnan(orientationOrder(i))
                    Screen('DrawTextures',w,gaborTexture{stimOrder(iBlock,i),2});
                    Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
                        4, [0,0,0], [rect(3)/2, rect(4)/2], 1);
                    changeTimesAbs = [changeTimesAbs time];
                end
            else % only fixation cross
                Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
                    4, [0,0,0], [rect(3)/2, rect(4)/2], 1);
            end
            
            T.presentationTime(iBlock,i,j) = Screen('Flip',w, time);
            time = time+1/flickerFreq;
            
            % Check for inputs
            keyPressed = false; % clear previous response
            while GetSecs < (time - 0.016) % look for responses until within one frame of a new presentation
                
                [~, keyTime, keyCode] = KbCheck(-3);
                key = KbName(keyCode);
                
                if ~iscell(key) % only start a keypress if there is only one key being pressed
                    if strcmp(key, resKey)
                        
                        % fill in B
                        respTime = keyTime - time;
                        R(iBlock,i,j) = respTime; % RT
                        
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
            
            if keyPressed
                responseCount = responseCount + 1;
                responseTimes(responseCount) = respTime;
            end
        end
        
        % If presTime is not a multiple of the flicker frequency, the
        % blocks will be slightly too short. Correct this.
        time = time + (presTime - floor(presTime*flickerFreq)/flickerFreq);
        
    end
    
    % End of block screen 
    if iBlock == round(nBlocks/2)
        totalBreakTime = round(30/TR)*TR; % make break duration a multiple of the TR.
                
        text = 'End of mini-block 1 out of 2! \n There will now be a 30s break. ';
        DrawFormattedText(w, text, 'center', 'center', [255 255 255]);
        Screen('Flip',w, time);
        time = time + 2;
        
        % 26 second break: empty screen
        Screen('Flip', w, time);
        time = time + (totalBreakTime - 4);
        
        % put the fixation dot back on the screen 2 seconds before the
        % break ends.
        Screen('DrawLines', w, [0 0 -10 10; -10 10 0 0],...
                4, [0,0,0], [rect(3)/2, rect(4)/2], 1);
        Screen('Flip', w, time);
        time = time + 2;
    end
    
end

save(fullfile(output,saveName)); % save everything

text = 'This is the end of this run, please wait while the scanner stops. ';
DrawFormattedText(w, text, 'center', 'center', [255 255 255]);
vbl = Screen('Flip', w);
WaitSecs(2);
Screen('CloseAll')
sca;
