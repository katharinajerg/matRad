function [cst, ct] = matRad_importDicomUSStructureSet(structureSetPath, imgPath, deformationDataPath)
% matRad_importDicomUSStructureSet is a matRad function to import 
% a predefined set of dicom files of an ultra sound with predefined contour
% data, which can be deformed or not, into matRad's native data formats
% 
% call
% [cst, ct] = matRad_importDicomUSStructureSet(structureSetPath, imgPath, deformationDataPath) 
%
% input
%   structureSetPath:    file to be imported which contains the structure set
%   imgPath:             file to be imported which contains one ultrasound image
%   deformationDataPath: file to be imported which contains deformation
%                        data (optional input)
%
% output
%   ct:     matRad ct struct
%   cst:    matRad cst struct
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 3
    performDeformation = 1;
elseif nargin == 2 
    performDeformation = 0;
end

info = dicominfo(structureSetPath);
imgInfo = dicominfo(imgPath); % US image contains information about the image geometry, which is not contained in the structure set
contour = dicomContours(info);
imgInfo.ImagePositionPatient(3) = contour.ROIs.ContourData{1,1}{end,1}(1,3); % set patient position to last slice, because the original US has more slices

figure
plotContour(contour)

% get geometric information
factor_x = 1/3; % number of pixel factor
factor_y = 1/3; % number of pixel factor
factor_z = 3;
imgSize = [floor(imgInfo.Rows*factor_y) floor(imgInfo.Columns*factor_x) factor_z*size(contour.ROIs.ContourData{1,1},1)]; % x = columns is second cube dimension
imgInfo.PixelSpacing(1) = imgInfo.PixelSpacing(1)/factor_x;
imgInfo.PixelSpacing(2) = imgInfo.PixelSpacing(2)/factor_y;
imgInfo.PixelSpacing(3) = (contour.ROIs.ContourData{1,1}{1,1}(1,3)-contour.ROIs.ContourData{1,1}{2,1}(1,3))/factor_z;

xlim = [imgInfo.ImagePositionPatient(1)-0.5*imgInfo.PixelSpacing(1), ...
        imgInfo.PixelSpacing(1)*double(imgSize(2)) + ...
        imgInfo.ImagePositionPatient(1)-0.5*imgInfo.PixelSpacing(1)];
ylim = [imgInfo.ImagePositionPatient(2)-0.5*imgInfo.PixelSpacing(2), ...
        imgInfo.PixelSpacing(2)*double(imgSize(1)) + ...
        imgInfo.ImagePositionPatient(2)-0.5*imgInfo.PixelSpacing(2)];
zlim = [contour.ROIs.ContourData{1,1}{end,1}(1,3), contour.ROIs.ContourData{1,1}{1,1}(1,3)];

referenceInfo = imref3d(imgSize,xlim,ylim,zlim);

roiName = contour.ROIs.Name;
cst = cell(height(roiName),6);
cst(:,2) = roiName;
completeMask = 0;
for i = 1:size(roiName)

    cst(i,1) = {i-1};
    if performDeformation == 1
        [mask, imgSize] = createDeformedMask(imgInfo, referenceInfo, ...
                    contour, imgSize, deformationDataPath, roiName(i));
    else
        mask = createMask(contour,char(roiName(i)),referenceInfo);
    end
    completeMask = completeMask + mask;
    cst{i,4}{1,1} = find(mask);
    cst{i,5} = struct('Priority',0,'Visible',1,'visibleColor', ...
        [0.33,0.667-1*0.2,1*0.2]);

end

% get geometric information
zPixelSpacing = imgInfo.PixelSpacing(3);

%% create ct struct

% use data from imported DICOM for most variables
dicomInfo = struct('PixelSpacing',imgInfo.PixelSpacing,'ImagePositionPatient' ...
    ,imgInfo.ImagePositionPatient,'ImageOrientationPatient', ...
    imgInfo.ImageOrientationPatient, 'PatientPosition', ...
    imgInfo.PatientPosition,'Width',imgInfo.Width,'Height', ...
    imgInfo.Height,'Manufacturer',imgInfo.Manufacturer,'ManufacturerModelName', ...
    imgInfo.ManufacturerModelName,'PatientName',imgInfo.PatientName);

% enter generated mask as ct cube
ct = struct('cube',{{double(completeMask)}},'resolution', ...
    struct('x',imgInfo.PixelSpacing(1),'y',imgInfo.PixelSpacing(2),'z', ...
    zPixelSpacing),'cubeDim',imgSize,'numOfCtScen',0,'dicomInfo', ...
    dicomInfo,'dicomMeta',info);

% determine x, y, z voxel coordinates of the mask in global coordinates
xLastVoxelPos = ct.dicomInfo.ImagePositionPatient(1) + (double(imgSize(2))-1) * imgInfo.PixelSpacing(1);
yLastVoxelPos = ct.dicomInfo.ImagePositionPatient(2) + (double(imgSize(1))-1) * imgInfo.PixelSpacing(2);
zLastVoxelPos = ct.dicomInfo.ImagePositionPatient(3) + (double(imgSize(3))-1) * imgInfo.PixelSpacing(3);

x = linspace(ct.dicomInfo.ImagePositionPatient(1),xLastVoxelPos,imgSize(2));
y = linspace(ct.dicomInfo.ImagePositionPatient(2),yLastVoxelPos,imgSize(1));
z = linspace(ct.dicomInfo.ImagePositionPatient(3),zLastVoxelPos,imgSize(3));

ct.x = x;
ct.y = y;
ct.z = z;

% for visualization the mask cube is used as replacement for the HU cube
ct.cubeHU = ct.cube;
end

%% function to create deformed masks

function [deformMask, imgSize] = createDeformedMask(imgInfo, referenceInfo, contour, imgSize, deformationDataPath, roiName)


originalMask = createMask(contour,char(roiName),referenceInfo);

% zero padding of 20 mm to account for deformed structure 
for i = 1:ceil(20/imgInfo.PixelSpacing(3)) 
    originalMask(:,:,size(originalMask,3)+1) = zeros([size(originalMask,1),size(originalMask,2)]);
end

tempSize = size(originalMask);
imgSize(3) = tempSize(3); % to account for zero padded slices
idx = find(originalMask);
[x,y,z] = ind2sub(size(originalMask), idx); % x is rows, y is colomns and z is slice

xLastVoxelPos = imgInfo.ImagePositionPatient(1) + (double(imgSize(2))-1) * imgInfo.PixelSpacing(1);
yLastVoxelPos = imgInfo.ImagePositionPatient(2) + (double(imgSize(1))-1) * imgInfo.PixelSpacing(2);
zLastVoxelPos = imgInfo.ImagePositionPatient(3) + (double(imgSize(3))-1) * imgInfo.PixelSpacing(3);

x_position = linspace(imgInfo.ImagePositionPatient(1),xLastVoxelPos,imgSize(2));
y_position = linspace(imgInfo.ImagePositionPatient(2),yLastVoxelPos,imgSize(1));
z_position = linspace(imgInfo.ImagePositionPatient(3),zLastVoxelPos,imgSize(3));

vox_pos_in_patient_space{1,1} = {};

for i = 1:numel(x)
    vox_pos_in_patient_space{i,1} = [x_position(x(i))];
    vox_pos_in_patient_space{i,2} = [y_position(y(i))];
    vox_pos_in_patient_space{i,3} = [z_position(z(i))];
end

% creating seprate variables for voxels in X, Y and Z for quicker
% calculation further in the code
voxX =  cell2mat(vox_pos_in_patient_space(:,1));
voxY =  cell2mat(vox_pos_in_patient_space(:,2));
voxZ =  cell2mat(vox_pos_in_patient_space(:,3));

%% Generate deformed contour
shifts = [-50,-70,-85];

% load deformation data and convert to mm
tissue(1) = struct(); 
imported_tissue = vtkRead(deformationDataPath);
tissue.verticies = 1000*imported_tissue.points;
tissue.deformation_field = 1000*imported_tissue.pointData.displacement;
clear imported_tissue

tissue.verticies(:,1) = tissue.verticies(:,1) + shifts(1); 
tissue.verticies(:,2) = tissue.verticies(:,2) + shifts(2);
tissue.verticies(:,3) = tissue.verticies(:,3) + shifts(3);

% remove verticies outside of physical domain
ind_zero_deformation = find(all(tissue.deformation_field < 1e-6,2));
ind_back_wall = find(tissue.verticies(:,3) > ((150+shifts(3))-1e-6));
ind_no_phase_field = ismember(ind_zero_deformation, ind_back_wall);
ind_phase_field = ind_zero_deformation;
ind_phase_field(ind_no_phase_field) = [];
tissue.verticies(ind_phase_field,:) = [];
tissue.deformation_field(ind_phase_field,:) = [];

% Remove doubled occurance of verticies.
[less_vert, ia, ic] = unique(tissue.verticies, 'stable', 'rows');
tissue.verticies = tissue.verticies(ia,:);
tissue.deformation_field = tissue.deformation_field(ia,:);

% interpolate deformation data with voxels in patient space
deformation_field(:,1) = griddata(double(tissue.verticies(:,1)), ...
    double(tissue.verticies(:,2)), double(tissue.verticies(:,3)),...
    double(tissue.deformation_field(:,1)), ...
    voxX, voxY, voxZ);
deformation_field(:,2) = griddata(double(tissue.verticies(:,1)), ...
    double(tissue.verticies(:,2)), double(tissue.verticies(:,3)),...
    double(tissue.deformation_field(:,2)), ...
    voxX, voxY, voxZ);
deformation_field(:,3) = griddata(double(tissue.verticies(:,1)), ...
    double(tissue.verticies(:,2)), double(tissue.verticies(:,3)),...
    double(tissue.deformation_field(:,3)), ...
    voxX, voxY, voxZ);

deformedField = deformation_field + cell2mat(vox_pos_in_patient_space(:,:));

%% find nearest neighbor mask points after deformation
deformMask = zeros(size(originalMask));
for i = 1:numel(x)
    deformedPoints = [deformedField(i,1), deformedField(i,2), deformedField(i,3)];
    diffInX = abs(x_position - deformedPoints(1));
    [diff_1,indX] = min(diffInX);
    diffInY = abs(y_position - deformedPoints(2));
    [diff_2,indY] = min(diffInY);
    diffInZ = abs(z_position - deformedPoints(3));
    [diff_3,indZ] = min(diffInZ);
    deformMask(indX,indY,indZ) = 1;
end

end
