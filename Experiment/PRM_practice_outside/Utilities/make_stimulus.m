function stimulus = make_stimulus(orientation,vis_level)

%%% --- Create the basic Gabors --- %%%
% rotation
rotAngle = -1 * (orientation+90);

% Gabor grating details
contrast = 1;
phase = 0;
spatialFrequency = 0.7;
gratingSizeDegrees = 5;
innerDegree = 0; %gratingSizeDegrees/15;

% Makes square gabor then masks with an outer and inner annulus to create a
% circular gabor with a hole for a fixation cross. Rotates the gabor to the
% desired angle.
[gaborPatch,~,annulusMatrix] = makeGabor(contrast, gratingSizeDegrees,...
    phase,spatialFrequency,innerDegree, rotAngle);

%%% --- Add noise to gabor --- %%%
if vis_level < 1
    stimulus = addNoise(gaborPatch, vis_level, annulusMatrix);
else
    stimulus = gaborPatch;
end
