function [nPixels, nPixelsUnrounded] = degrees2pixels(degrees)
% [nPixels, nPixelsUnrounded] = degrees2pixels(degrees, distFromScreen_inCm, pixels_perCm, screenNumber)
%
% Converts degrees to pixels, given distance from screen in centimeters
% and pixels per centimeter. Result is rounded by default. Use
% nPixelsUnrounded to get the unrounded result of the calculation.
%

% global distFromScreen pixelsPerCm;
distFromScreen = 60; 
pixelsPerCm = (1920/56 + 1200/36.5)/2;

% convert degrees to centimeters
sizeOnScreen_inCm = 2 * distFromScreen * tan((degrees/2) * (pi/180));

% convert centimeters to pixels
nPixelsUnrounded = sizeOnScreen_inCm * pixelsPerCm;
nPixels = round(nPixelsUnrounded);