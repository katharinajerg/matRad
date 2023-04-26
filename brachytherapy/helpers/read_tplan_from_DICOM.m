clear all
close all

% Pat
for patient = 56:117
%     path = ['~/thindrives/ProstateData/',num2str(patient),'/IntraOp/IntraOp/'];
    path = ['~/thindrives/ProstateData/Pat',num2str(patient),'/'];
    D = '~/thindrives/2023_01_09_LDR_tissue_deformation/dicom-dict-iotp.txt';
    info_pl = dicominfo([path,'PL001.dcm'],'dictionary',D);
    
    % source locations needle insertion
    seeds = struct2array(info_pl.ApplicationSetupSequence);
    tplan = cell(length(seeds),3);
    
    for i = 1:length(seeds)
        tplan{i,1} = seeds(i).ApplicationSetupNumber;
        tplan{i,2} = seeds(i).ChannelSequence.Item_1.BrachyControlPointSequence.Item_1.ControlPoint3DPosition;
        tplan{i,3} = seeds(i).ChannelSequence.Item_1.BrachyControlPointSequence.Item_1.ControlPointOrientation;
    end
    
    
    save([path,'tplan_orig.mat'], "tplan")


end
