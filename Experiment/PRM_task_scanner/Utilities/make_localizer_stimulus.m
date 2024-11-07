function stimulus = make_localizer_stimulus(orientation,colour)

%%% --- Create the basic Gabors --- %%%
% rotation
rotAngle = -1 * (orientation+90);

% Gabor grating details
contrast = 1;
phase    = 0;
spatialFrequency = 0.7;
gratingSizeDegrees = 5;
innerDegree = 0; %gratingSizeDegrees/15;

% Makes square gabor then masks with an outer and inner annulus to create a
% circular gabor with a hole for a fixation cross. Rotates the gabor to the
% desired angle.
[gaborPatch,~,annulusMatrix] = makeGabor(contrast, gratingSizeDegrees,...
    phase,spatialFrequency,innerDegree, rotAngle);

if colour > 0 
    tmpStim = repmat(gaborPatch,1,1,3); % to create rgb channels
    rgb = 1:3; 
    for r = 1:length(rgb)
        if r~=colour % turn off colour stimulus for non-coloured one
            tmp = gaborPatch; tmp(annulusMatrix>0) = tmp(annulusMatrix>0).*0.7;
            tmpStim(:,:,r) = tmp;
        end
    end
    stimulus = im2double(tmpStim);
else
stimulus = gaborPatch;
end