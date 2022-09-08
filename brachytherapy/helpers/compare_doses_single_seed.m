clear all;
close all;

% params
slice_num = 4;
percent = 1;
dta = 1; % mm
local = 0; % Perform global gamma

% load matRad dose
load("result_100.mat")
matrad_dose_100 = resultGUI.physicalDose;
load("result_010.mat")
matrad_dose_010 = resultGUI.physicalDose;
load("dij.mat")
matrad_cube = dij.ctGrid;

% calc difference
diff = matrad_dose_010 - matrad_dose_100;

%% plot
figure
imagesc(diff(:,:,4))
title('diff')
colorbar
figure
imagesc(matrad_dose_100(:,:,slice_num))
title('dose [100]')
colorbar
figure
imagesc(matrad_dose_010(:,:,slice_num))
title('dose [010]')
colorbar

