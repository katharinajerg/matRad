%read DICOM file
clear all;
close all;
clc;

patient = 90;
path_1 = ['~/thindrives/ProstateData/',num2str(patient),'/IntraOp/IntraOp/'];

D = '/home/kjerg/Daten/dicom-dict-iotp.txt';
info_pl_1 = dicominfo([path_1,'PL001.dcm'],'dictionary',D);
info_ss_1 = dicominfo([path_1,'SS001.dcm']);
pathStructureSet = [path_1, 'SS001.dcm']; 
pathImg = [path_1, 'MR001.dcm'];


% source locations plan
seeds_1 = struct2array(info_pl_1.ApplicationSetupSequence);
tplan_1 = cell(length(seeds_1),3);

for i = 1:length(seeds_1)
    tplan_1{i,1} = seeds_1(i).ApplicationSetupNumber;
    tplan_1{i,2} = seeds_1(i).ChannelSequence.Item_1.BrachyControlPointSequence.Item_1.ControlPoint3DPosition;
    tplan_1{i,3} = seeds_1(i).ChannelSequence.Item_1.BrachyControlPointSequence.Item_1.ControlPointOrientation;
end


% % strucure set
contour_1 = dicomContours(info_ss_1);

% cont = struct(contour_1);
% cont.ROIs.Color{1}(2) = 200;
% cont.ROIs.Color{2}(2) = 255;
% cont.ROIs.Color{3}(2) = 200;
% cont.ROIs.Number(1) = 3;
% out = addContour(contour_1, cont.ROIs.Number(1), cont.ROIs.Name{1}, cont.ROIs.ContourData{1},cont.ROIs.GeometricType{1},cont.ROIs.Color{1});   
% out2 = addContour(out, 4, cont.ROIs.Name{2}, cont.ROIs.ContourData{2},cont.ROIs.GeometricType{2},cont.ROIs.Color{2});   
% out3 = addContour(out2, 5, cont.ROIs.Name{3}, cont.ROIs.ContourData{3},cont.ROIs.GeometricType{3},cont.ROIs.Color{3});   


% create surface
[cst, ct] = matRad_importDicomUSStructureSet(pathStructureSet,pathImg);

% Prostate bed objective
cst{1,3} = 'TARGET';
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(400,200));
cst{1,5}.Priority = 3;

% find all target voxels from cst cell array
V = [];
for i=1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        V = [V;vertcat(cst{i,4}{:})];
    end
end

% Remove double voxels
V = unique(V);
% generate voi cube for targets
voiTarget    = zeros(ct.cubeDim);
voiTarget(V) = 1;
    
% throw error message if no target is found
if isempty(V)
    matRad_cfg.dispError('Could not find target.');
end

% Convert linear indices to 3D voxel coordinates
[coordsY_vox, coordsX_vox, coordsZ_vox] = ind2sub(ct.cubeDim,V);


% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));

% ignore densities outside of contours
eraseCtDensMask = ones(prod(ct.cubeDim),1);
eraseCtDensMask(V) = 0;
for i = 1:ct.numOfCtScen
    ct.cube{i}(eraseCtDensMask == 1) = 0;
end

TargX = ct.x(coordsX_vox);
TargY = ct.y(coordsY_vox);
TargZ = ct.z(coordsZ_vox);
%Prostate = plot3(TargX,TargY,TargZ,'.', 'Color','b','DisplayName', 'prostate');

P = [TargX',TargY',TargZ'];
k = boundary(P,1);


% plot
figure(1)
hold on 
for i = 1:length(tplan_1)
    p = scatter3(tplan_1{i,2}(1), tplan_1{i,2}(2), tplan_1 {i,2}(3), 'g', 'filled');
end
plotContour(contour_1);
% plotContour(out3);

trisurf(k,P(:,1),P(:,2),P(:,3),'FaceColor','red','FaceAlpha',0.1,'LineStyle','none')