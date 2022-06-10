function [cst, ct] = matRad_importDicomUSStructureSet(structureSetPath, imgPath)
% matRad_importDicomUSStructureSet is a matRad function to import 
% a predefined set of dicom files of an ultra sound with predefined contour
% data into matRad's native data formats
% 
% call
% [cst, ct] = matRad_importDicomUSStructureSet(structureSetPath, imgPath)  
% [cst, ct] = matRad_importDicomUSStructureSet(structureSetPath, imgPath)
%
% input
%   structureSetPath: file to be imported which contains the structure set
%   imgPath:          file to be imported which contains one ultrasound image
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
arguments 
    structureSetPath {mustBeFile}
    imgPath {mustBeFile}
end

%% read DICOM data
info = dicominfo(structureSetPath);
imgInfo = dicominfo(imgPath); % US image contains information about the image geometry, which is not contained in the structure set
contour = dicomContours(info);
imgInfo.ImagePositionPatient(3) = contour.ROIs.ContourData{1,1}{end,1}(1,3); % set patient position to last slice, because the original US has more slices

figure
plotContour(contour)

%% create masks from the contour data and convert to linear indices
% get geometric information
imgSize = [imgInfo.Rows imgInfo.Columns size(contour.ROIs.ContourData{1,1},1)]; % x = columns is second cube dimension
zPixelSpacing = contour.ROIs.ContourData{1,1}{1,1}(1,3)-contour.ROIs.ContourData{1,1}{2,1}(1,3);

xlim = [imgInfo.ImagePositionPatient(1)-0.5*imgInfo.PixelSpacing(1), ...
    imgInfo.PixelSpacing(1)*double(imgInfo.Columns) + imgInfo.ImagePositionPatient(1)-0.5*imgInfo.PixelSpacing(1)];
ylim = [imgInfo.ImagePositionPatient(2)-0.5*imgInfo.PixelSpacing(2), ...
    imgInfo.PixelSpacing(2)*double(imgInfo.Rows) + imgInfo.ImagePositionPatient(2)-0.5*imgInfo.PixelSpacing(2)];
zlim = [contour.ROIs.ContourData{1,1}{end,1}(1,3), contour.ROIs.ContourData{1,1}{1,1}(1,3)];

referenceInfo = imref3d(imgSize,xlim,ylim,zlim);

% for all ROIs create a mask and convert to linear indicies
roiName = contour.ROIs.Name;
cst = cell(height(roiName),6);
cst(:,2) = roiName;
mk = 0;
for i = 1:size(roiName)

    cst(i,1) = {i-1};
    mask = createMask(contour,char(roiName(i)),referenceInfo);
    mk = mk + mask;
    % generate the indices using the 'find' function which has in-built 
    % ind2sub function
    temp = find(mask);
    cst{i,4}{1,1} = temp;
    cst{i,5} = struct('Priority',0,'Visible',1,'visibleColor', ...
        [0.33,0.667-i*0.2,i*0.2]);

end

%% create ct struct
% use data from imported DICOM for most variables
dicomInfo = struct('PixelSpacing',imgInfo.PixelSpacing,'ImagePositionPatient' ...
    ,imgInfo.ImagePositionPatient,'ImageOrientationPatient', ...
    imgInfo.ImageOrientationPatient, 'PatientPosition', ...
    imgInfo.PatientPosition,'Width',imgInfo.Width,'Height', ...
    imgInfo.Height,'Manufacturer',imgInfo.Manufacturer,'ManufacturerModelName', ...
    imgInfo.ManufacturerModelName,'PatientName',imgInfo.PatientName);

% enter generated mask as ct cube
ct = struct('cube',{{double(mk)}},'resolution', ...
    struct('x',imgInfo.PixelSpacing(1),'y',imgInfo.PixelSpacing(2),'z', ...
    zPixelSpacing),'cubeDim',imgSize,'numOfCtScen',0,'dicomInfo', ...
    dicomInfo,'dicomMeta',info);

% determine x, y, z voxel coordinates of the mask in global coordinates
xLastVoxelPos = ct.dicomInfo.ImagePositionPatient(1) + (double(imgSize(2))-1) * ct.dicomInfo.PixelSpacing(1);
yLastVoxelPos = ct.dicomInfo.ImagePositionPatient(2) + (double(imgSize(1))-1) * ct.dicomInfo.PixelSpacing(2);
zLastVoxelPos = ct.dicomInfo.ImagePositionPatient(3) + (double(imgSize(3))-1) * zPixelSpacing;

x = linspace(ct.dicomInfo.ImagePositionPatient(1),xLastVoxelPos,imgSize(2));
y = linspace(ct.dicomInfo.ImagePositionPatient(2),yLastVoxelPos,imgSize(1));
z = linspace(ct.dicomInfo.ImagePositionPatient(3),zLastVoxelPos,imgSize(3));

ct.x = x;
ct.y = y;
ct.z = z;

% for visualization the mask cube is used as replacement for the HU cube
ct.cubeHU = ct.cube;
end
