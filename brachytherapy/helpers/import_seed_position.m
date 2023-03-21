%read DICOM file
% clear all;
% close all;
% clc;


path_1 = '~/Daten/Pat1/needle-insertion/';

D = '/home/kjerg/Daten/dicom-dict-iotp.txt';
info_pl_1 = dicominfo([path_1,'PL001.dcm'],'dictionary',D);
info_ss_1 = dicominfo([path_1,'SS001.dcm']);


% source locations plan
seeds_1 = struct2array(info_pl_1.ApplicationSetupSequence);
tplan_1 = cell(length(seeds_1),3);

for i = 1:length(seeds_1)
    tplan_1{i,1} = seeds_1(i).ApplicationSetupNumber;
    tplan_1{i,2} = seeds_1(i).ChannelSequence.Item_1.BrachyControlPointSequence.Item_1.ControlPoint3DPosition;
    tplan_1{i,3} = seeds_1(i).ChannelSequence.Item_1.BrachyControlPointSequence.Item_1.ControlPointOrientation;
end

figure
hold on 
for i = 1:length(tplan_1)
    p = scatter3(tplan_1{i,2}(1), tplan_1{i,2}(2), tplan_1 {i,2}(3), 'g');
end


% % strucure set
contour_1 = dicomContours(info_ss_1);

cont = struct(contour_1);
cont.ROIs.Color{1}(2) = 200;
cont.ROIs.Color{2}(2) = 255;
cont.ROIs.Color{3}(2) = 200;
cont.ROIs.Number(1) = 3;
out = addContour(contour_1, cont.ROIs.Number(1), cont.ROIs.Name{1}, cont.ROIs.ContourData{1},cont.ROIs.GeometricType{1},cont.ROIs.Color{1});   
out2 = addContour(out, 4, cont.ROIs.Name{2}, cont.ROIs.ContourData{2},cont.ROIs.GeometricType{2},cont.ROIs.Color{2});   
out3 = addContour(out2, 5, cont.ROIs.Name{3}, cont.ROIs.ContourData{3},cont.ROIs.GeometricType{3},cont.ROIs.Color{3});   

hold on 
plotContour(out3);
