function [gaborPatch, widthArray,annulusMatrix] = makeGabor(contrastGrating,sizeGrating,phaseGrating,spatialFrequency, sizeInnerCircleDegrees, angle)
% The grating is first created in terms of desired luminance values, and
% then these are converted to RGB. Slow, but (relatively) clean.
% spatialFrequency: cycles/degree

%tic

width = degrees2pixels(sizeGrating);

sizeInnerCircle = degrees2pixels(sizeInnerCircleDegrees);
%startLinearDecayDegrees = 1;
startLinearDecayDegrees = sizeGrating/15;
startLinearDecay = degrees2pixels(startLinearDecayDegrees);

nCycles = sizeGrating*spatialFrequency; % number of cycles in a stimulus

% compute the pixels per grating period
pixelsPerGratingPeriod = width / nCycles;

spatialFrequency = 1 / pixelsPerGratingPeriod; % How many periods/cycles are there in a pixel?
radiansPerPixel = spatialFrequency * (2 * pi); % = (periods per pixel) * (2 pi radians per period)

% Make luminance calculations based on the grating contrast value (this will change during staircasing)
[~, ~, ~, G_Lmin, G_Lmax, lum] = calibrateLum(contrastGrating);
G_background = (G_Lmin+G_Lmax)/2;
G_lumRange = G_Lmax - G_background;

halfWidthOfGrid = width / 2;
widthArray = (-halfWidthOfGrid) : halfWidthOfGrid;  % widthArray is used in creating the meshgrid.

% Creates a two-dimensional square grid.  For each element i = i(x0, y0) of
% the grid, x = x(x0, y0) corresponds to the x-coordinate of element "i"
% and y = y(x0, y0) corresponds to the y-coordinate of element "i"
[x y] = meshgrid(widthArray);

% make annulusMatrix (without using for loops)
stimRadii      = sqrt(x.^2 + y.^2);
annulusMatrix = stimRadii <= (width+1)/2;

% Creates a sinusoidal grating, where the period of the sinusoid is
% approximately equal to "pixelsPerGratingPeriod" pixels.
% Note that each entry of gratingMatrix varies between minus one and
% one; -1 <= gratingMatrix(x0, y0)  <= 1
stimulusMatrix = sin(radiansPerPixel * x + phaseGrating);
stimulusMatrix = imrotate(stimulusMatrix, angle); 

% Make a fading annulus, to use as a mask.
annulusMatrix = makeLinearMaskCircleAnn(width+1,width+1,sizeInnerCircle,startLinearDecay,width/2+1);
annulusMatrix = imrotate(annulusMatrix, 45);
sSm = size(stimulusMatrix,1); sAm = size(annulusMatrix,1); sDf = sAm-sSm;
if sDf > 0 % to make sure annulus and stim have the same size
    stimTmp        = zeros(size(annulusMatrix));
    stimTmp(sDf/2+1:sAm-(sDf/2),sDf/2+1:sAm-(sDf/2)) = stimulusMatrix;
    stimulusMatrix = stimTmp;
end
stimulusMatrix = stimulusMatrix .* annulusMatrix;

% create luminance-defined grating
gaborPatch = G_lumRange*stimulusMatrix + G_background;

% the gaborPatch is currently defined in terms of luminance; we need to
% convert it to RGB.
gaborVect = reshape(gaborPatch,numel(gaborPatch),1);
gaborVect = interp1(lum,0:255,gaborVect,'nearest');
gaborPatch = reshape(gaborVect,size(gaborPatch));

end
