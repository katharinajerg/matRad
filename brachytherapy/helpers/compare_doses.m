clear all;
close all;

% params
slice_num = 5;
percent = 3;
dta = 3; % mm
local = 0; % Perform global gamma

% read DICOM infos
info_do = dicominfo('~/Daten/DO001.dcm');
info_pl = dicominfo('~/Daten/PL001.dcm');
info_mr = dicominfo('~/Daten/MR001.dcm');

% read dose
do = dicomread(info_do);
dose = squeeze(do);
% delete first two and last two rows, because they are not imported and
% therefore calculated in matRad
dose = dose(:,:,3:end-2);
% scale to GY
dose = info_do.DoseGridScaling*double(dose);

% load matRad dose
load("result.mat")
matrad_dose = resultGUI.physicalDose;
load("dij.mat")
matrad_cube = dij.ctGrid;

% transform intrinsic dose points from voxels to patient coordinates 
pixel_spacing_do = info_do.PixelSpacing;
image_position_do = info_do.ImagePositionPatient;
image_orientation_do = info_do.ImageOrientationPatient;

dose_cube.x = image_position_do(1):pixel_spacing_do:image_position_do(1)+(size(dose,2)-1)*pixel_spacing_do(1);
dose_cube.y = image_position_do(2):pixel_spacing_do:image_position_do(2)+(size(dose,1)-1)*pixel_spacing_do(2);
dose_cube.z = flip(dij.ctGrid.z);

% interpolate fine matrad grid to coarse variseed grid
interpolated_dose = interp3(matrad_cube.x, matrad_cube.y', matrad_cube.z, ...
    matrad_dose, dose_cube.x,dose_cube.y',dose_cube.z,'linear',0);

% cutoff large values
interpolated_dose(interpolated_dose > 700) = 700;

% cutoff boundary which doesn't contain dose values in MatRad
[rows, cols, values] = find(interpolated_dose(:,:,1));
interpolated_dose = interpolated_dose(min(rows):max(rows), min(cols):max(cols),:);
dose = dose(min(rows):max(rows), min(cols):max(cols),:);
dose_cube.x = dose_cube.x(min(rows):max(rows));
dose_cube.y = dose_cube.y(min(cols):max(cols));

% prepare data for gamma index calculation
reference.start = [dose_cube.x(1) dose_cube.y(1)]; % mm
reference.width = [dose_cube.x(2)-dose_cube.x(1) dose_cube.y(2)-dose_cube.y(1)]; % mm
reference.data = dose(:,:,slice_num);

target.start = [dose_cube.x(1) dose_cube.y(1)]; % mm
target.width = [dose_cube.x(2)-dose_cube.x(1) dose_cube.y(2)-dose_cube.y(1)]; % mm
target.data = interpolated_dose(:,:,slice_num);

% calc gamma index
gamma = CalcGamma(reference, target, percent, dta, 'local', local);
mean_gamma = mean(gamma(:))


%% plot doses
% figure(1)
% sliceViewer(dose, 'DisplayRange', [0,300]);
% figure(2)
% sliceViewer(interpolated_dose, 'DisplayRange', [0 300]);
% figure(3)
% sliceViewer(diff, 'DisplayRange', [0 50]);
% mean(diff(:))
% max(diff(:))
% min(diff(:))
figure
imagesc(dose(:,:,slice_num))
title('VariSeed dose')
colorbar
figure
imagesc(interpolated_dose(:,:,slice_num))
title('MatRad dose')
colorbar
figure
imagesc(gamma)
title(['gamma index (',num2str(percent),'%,',num2str(dta),'mm)'])
colorbar

