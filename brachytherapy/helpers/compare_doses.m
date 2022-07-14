clear all;
close all;

% params
slice_num = 5;
percent = 3;
dta = 3; % mm
local = 0; % Perform global gamma

%% read DICOM infos
info_do = dicominfo('~/Daten/Pat4/DO001.dcm');
info_pl = dicominfo('~/Daten/Pat4/PL001.dcm');
info_mr = dicominfo('~/Daten/Pat4/MR001.dcm');

%% read dose
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

%% get dosimetry measures from imported plan
num_points = size(info_do.DVHSequence.Item_1.DVHData,1)/2;
bin_size = info_do.DVHSequence.Item_1.DVHData(1);
for i = 1:num_points
    v1(i) = info_do.DVHSequence.Item_1.DVHData(2*i);
    v2(i) = info_do.DVHSequence.Item_2.DVHData(2*i);
    v3(i) = info_do.DVHSequence.Item_3.DVHData(2*i);
end
d = bin_size:bin_size:num_points*bin_size;
v_rel_1 = v1/v1(1)*100;
v_rel_2 = v2/v2(1)*100;
v_rel_3 = v3/v3(1)*100;

figure
hold on 
plot(d,v_rel_1)
plot(d,v_rel_2)
plot(d,v_rel_3)
xlabel('dose in Gy')
ylabel('relative volume in %')
target_dose = 160;
P_D90  = getDx(v_rel_1, 90, bin_size);
P_V100 = v_rel_1(target_dose/bin_size);
P_V150 = v_rel_1(1.5*target_dose/bin_size);
U_D10  = getDx(v_rel_2, 10, bin_size);
U_D30  = getDx(v_rel_2, 30, bin_size);
R_D2cc = getDx(v3, 2, bin_size);
R_D01cc = getDx(v3, 0.1, bin_size);
R_V100 = v_rel_3(target_dose/bin_size);

fprintf('\n\tdosimetric parameters\n')
fprintf(1 + ~(P_D90 > 0.9*target_dose), 'P_D90 is %3.2f Gy and should be greater than %d Gy.\n', P_D90, 0.9*target_dose);
fprintf(1 + ~(P_V100 > 95), 'P_V100 is %3.2f%% and should be greater than %d%%.\n', P_V100, 95);
fprintf(1 + ~(P_V150 > 45 && P_V150 < 65), 'P_V150 is %3.2f%% and should be between %d%% and %d%%.\n', P_V150, 45, 65);
fprintf(1 + ~(U_D10 < 1.5*target_dose), 'U_D10 is %3.2f Gy and should be less than %d Gy.\n', U_D10, 1.5*target_dose);
fprintf(1 + ~(U_D30 < 1.3*target_dose), 'P_D30 is %3.2f Gy and should be less than %d Gy.\n', U_D30, 1.3*target_dose);
fprintf(1 + ~(R_D2cc < 145), 'R_D2cc is %3.2f Gy and should be less than %d Gy.\n', R_D2cc, 145);
fprintf(1 + ~(R_D01cc < 200), 'R_D01cc is %3.2f Gy and should be less than %d Gy.\n\n', R_D01cc, 200);

function Dx = getDx(v_rel, x, bin_size)
    for i = 1:numel(v_rel)
        if (v_rel(i) < x)
            Dx = (i-1)*bin_size + (v_rel(i-1)-x)/(v_rel(i-1)-v_rel(i))*bin_size;
            break;
        end
    end
end

