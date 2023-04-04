function [cst, ct] = matRad_importDicomDummySphericalStructureSet(structureSetPath, imgPath)
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

contour = contour.deleteContour(1);
contour = contour.deleteContour(2);
surfaces = cell(6,1);
surface = zeros(16,3);%cell2mat(contour.ROIs(1,:).ContourData{1,1}(1,1));
r = 3;
m = [0,0,0];
z = linspace(-r+0.1,r-0.1,10);
z=flip(z);
for i = 1:numel(z)
    surface(:,1) = 0;
    surface(:,2) = 0;
    surface(:,3) = z(i);
    
    y12 = sqrt(r^2 - z(i)^2);
    xb = 0.5*sqrt(y12^2/2);
    yb = sqrt(r^2 -z(i)^2 - 1/8 * y12^2);

    surface(1,1) = 0;
    surface(1,2) = y12;
    surface(2,1) = -xb;
    surface(2,2) = yb;
    surface(3,1) = -sqrt(y12^2/2);
    surface(3,2) = sqrt(y12^2/2);
    surface(4,1) = -yb;
    surface(4,2) = xb;
    surface(5,1) = -y12;
    surface(5,2) = 0;   
    surface(6,1) = -yb;
    surface(6,2) = -xb;
    surface(7,1) = -sqrt(y12^2/2);
    surface(7,2) = -sqrt(y12^2/2);
    surface(8,1) = -xb;
    surface(8,2) = -yb;
    surface(9,1) = 0;
    surface(9,2) = -y12;
    surface(10,1) = xb;
    surface(10,2) = -yb;
    surface(11,1) = sqrt(y12^2/2);
    surface(11,2) = -sqrt(y12^2/2);
    surface(12,1) = yb;
    surface(12,2) = -xb;
    surface(13,1) = y12;
    surface(13,2) = 0;
    surface(14,1) = yb;
    surface(14,2) = xb;
    surface(15,1) = sqrt(y12^2/2);
    surface(15,2) = sqrt(y12^2/2);
    surface(16,1) = xb;
    surface(16,2) = yb;

    surface(:,1) = surface(:,1) + m(1);
    surface(:,2) = surface(:,2) + m(2);
    surface(:,3) = surface(:,3) + m(3);

    surfaces(i,1) = mat2cell(surface,[16]);
end



geometricType = contour.ROIs.GeometricType{1};
for i = 1:numel(z)
    geometricType(i,1) = geometricType(1,1);
end

contour = addContour(contour, 2, contour.ROIs.Name{1}, surfaces, geometricType, contour.ROIs.Color{1});   

contour = contour.deleteContour(0);

imgInfo.ImagePositionPatient(3) = contour.ROIs.ContourData{1,1}{end,1}(1,3); % set patient position to last slice, because the original US has more slices

figure
plotContour(contour)

%% create masks from the contour data and convert to linear indices
% get geometric information
xlim = [-r, r];
ylim = [-r, r];
zlim = [-r, r];
imgSize = repmat(20, [1,3]);

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
    ,[xlim(1);ylim(1);zlim(1)],'ImageOrientationPatient', ...
    imgInfo.ImageOrientationPatient, 'PatientPosition', ...
    imgInfo.PatientPosition,'Width',imgInfo.Width,'Height', ...
    imgInfo.Height,'Manufacturer',imgInfo.Manufacturer,'ManufacturerModelName', ...
    imgInfo.ManufacturerModelName,'PatientName',imgInfo.PatientName);

% determine x, y, z voxel coordinates of the mask in global coordinates
x = linspace(xlim(1),xlim(2),imgSize(2));
y = linspace(ylim(1),ylim(2),imgSize(1));
z = linspace(zlim(1),zlim(2),imgSize(3));

% enter generated mask as ct cube
ct = struct('cube',{{double(mk)}},'resolution', ...
    struct('x',x(2)-x(1),'y',y(2)-y(1),'z', ...
    z(2)-z(1)),'cubeDim',imgSize,'numOfCtScen',0,'dicomInfo', ...
    dicomInfo,'dicomMeta',info);


ct.x = x;
ct.y = y;
ct.z = z;

% for visualization the mask cube is used as replacement for the HU cube
ct.cubeHU = ct.cube;
end
