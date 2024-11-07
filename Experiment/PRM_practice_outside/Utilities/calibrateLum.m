% This version of calibrate luminance is adjusted to include BOTH rgb and
% luminance terms in the calibration of a given contrast value - this
% calculates the luminance for every rgb value.
function [background, Lmin_rgb, Lmax_rgb, Lmin, Lmax, lum] = calibrateLum(contrast)

global environment scanner;

rgb =          [0        30      60      90      120     150     180     190     210     220     230     240     250     255];

if strcmp(environment,'mri') || strcmp(environment,'mri_offline')
    % made-up luminance values for the MRI projector.
        rgb         = [0        30      60      90      120     150     180     190     200     210     220     230     240     250     255];
        lum(1,:)    = [0.7      13.6    78.5    192     333     490     667     732     784     833     890     960     1000    1045    1060];
    lum = mean(lum,1);
elseif strcmp(environment,'workstation')
    % use made-up values for my workstation.
    lum(1,:) =     [0.27  4.1    24.6   63.8   109     163     224     243     286     302     317     325     326     327];
	lum = mean(lum,1);
elseif strcmp(environment,'macbook')
	% use made-up values for my macbook.
	lum(1,:) =     [0.27  4.1    24.6   63.8   109     163     224     243     286     302     317     325     326     327];
	lum = mean(lum,1);
else
    % not a calibrated environment; use beh1 lab as a proxy.
    disp('Warning: no luminance values found!');
    lum(1,:) =     rgb;%[0.2      2.30    9       18      35      53      80      91      110     115     132     136     151     156];
    lum = mean(lum,1);
end

% interpolate over the whole rgb range
lum = interp1(rgb,lum,0:255,'pchip');

%middle of the luminance range:
medium_lum = (min(lum) + max(lum))/2;
lum_diff = abs(lum - medium_lum);
[val ind] = min(lum_diff);
background = ind - 1; %ind is 1-256 and need 0-255

if exist('contrast', 'var')
    Lmin = medium_lum - (medium_lum-min(lum))*contrast;
    Lmax = medium_lum + (max(lum)-medium_lum)*contrast;
    
    lum_diff = abs(lum - Lmin);
    [val ind] = min(lum_diff);
    Lmin_rgb = ind - 1; %ind is 1-256 and need 0-255
    
    lum_diff = abs(lum - Lmax);
    [val ind] = min(lum_diff);
    Lmax_rgb = ind - 1; %ind is 1-256 and need 0-255    
end

%plot(0:255,lum)
