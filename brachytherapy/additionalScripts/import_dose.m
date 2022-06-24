%read DICOM file
clear all;
close all;

% read DICOM infos
info_do = dicominfo('~/Daten/DO001.dcm');
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

%% interpolate fine matrad grid to coarse variseed grid
interpolated_dose = interp3(matrad_cube.x, matrad_cube.y', matrad_cube.z, ...
    matrad_dose, dose_cube.x,dose_cube.y',dose_cube.z,'linear',0);

% cutoff large values
interpolated_dose(interpolated_dose > 700) = 700;

%% calculate difference
diff = double(dose) - interpolated_dose;


%% plot doses
figure(1)
sliceViewer(dose, 'DisplayRange', [0,300]);
figure(2)
sliceViewer(interpolated_dose, 'DisplayRange', [0 300]);
figure(3)
sliceViewer(diff, 'DisplayRange', [0 50]);
mean(diff(:))
max(diff(:))
min(diff(:))
figure
imagesc(dose(:,:,4))
colorbar
figure
imagesc(interpolated_dose(:,:,4))
colorbar
figure
imagesc(diff(:,:,4))
colorbar

