function [cst, ct] = matRad_importDicomDummyStructureSet(structureSetPath, imgPath)
% matRad_importDicomDummyStructureSet is a matRad function to create 
% a predefined set of dicom files of an ultra sound with rectangular contour
% data into matRad's native data formats
% 
% call
% [cst, ct] = matRad_importDicomDummyStructureSet(structureSetPath, imgPath)  
% [cst, ct] = matRad_importDicomDummyStructureSet(structureSetPath, imgPath)
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
testTarget = [-15, -45, -10; 10, -45, -10; 10, -20, -10; -15, -20, -10; -15,-45,-10;
              -15, -45, -15; 10, -45, -15; 10, -20, -15; -15, -20, -15; -15,-45,-15;
              -15, -45, -20; 10, -45, -20; 10, -20, -20; -15, -20, -20; -15,-45,-20; 
              -15, -45, -25; 10, -45, -25; 10, -20, -25; -15, -20, -25; -15,-45,-25;
              -15, -45, -30; 10, -45, -30; 10, -20, -30; -15, -20, -30; -15,-45,-30;];
testTarget = mat2cell(testTarget, [5 5 5 5 5]);

testOAR1     = [10, -45, -10; 20, -45, -10; 20, -20, -10; 10, -20, -10; 10,-45,-10;
                10, -45, -15; 20, -45, -15; 20, -20, -15; 10, -20, -15; 10,-45,-15;
                10, -45, -20; 20, -45, -20; 20, -20, -20; 10, -20, -20; 10,-45,-20; 
                10, -45, -25; 20, -45, -25; 20, -20, -25; 10, -20, -25; 10,-45,-25;
                10, -45, -30; 20, -45, -30; 20, -20, -30; 10, -20, -30; 10,-45,-30;];
testOAR1 = mat2cell(testOAR1, [5 5 5 5 5]);

testOAR2     = [-25, -45, -10; -15, -45, -10; -15, -20, -10; -25, -20, -10; -25, -45, -10;
                -25, -45, -15; -15, -45, -15; -15, -20, -15; -25, -20, -15; -25, -45, -15;
                -25, -45, -20; -15, -45, -20; -15, -20, -20; -25, -20, -20; -25, -45, -20; 
                -25, -45, -25; -15, -45, -25; -15, -20, -25; -25, -20, -25; -25, -45, -25; 
                -25, -45, -30; -15, -45, -30; -15, -20, -30; -25, -20, -30; -25, -45, -30;];
testOAR2= mat2cell(testOAR2, [5 5 5 5 5]);

contour= deleteContour(contour,0);
contour= deleteContour(contour,1);
contour= deleteContour(contour,2);
contour= addContour(contour,0,'TestTarget',testTarget,'CLOSED_PLANAR',[255;0;0]);
contour= addContour(contour,1,'TestOAR1',testOAR1,'CLOSED_PLANAR',[0;150;0]);
contour= addContour(contour,2,'TestOAR2',testOAR2,'CLOSED_PLANAR',[0;150;0]);

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
