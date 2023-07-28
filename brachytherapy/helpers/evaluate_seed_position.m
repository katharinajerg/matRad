% 14.03.2023, Katharina Jerg
% This script is used to deterimine the impact of each seed on different
% dose measures. 

clear all
close all

% MAKE SURE TO DEFINE:
% 1. in brachytherapy\matRad_automaticDifferentiation.m, l.62ff: wanted dose parameter
% 2. l.39f, folder for each dose parameter

% 3. import seed positions corresponding to the patient number
for patient = 1:35
    patientId = 1000+patient;
    load(['..\BRACHYTHERAPY_data\',num2str(patient),'\IntraOp\IntraOp\tplan_orig.mat']);
    numSeeds = size(tplan, 1);
    
    allSeeds = [];
    for i = 1:numSeeds
        if (false)%any(i == maxSeedIds)) %change position of one single seed
            pos = cell2mat(tplan(i,2))';
            pos(2) = pos(2)-20;
            allSeeds = [allSeeds,pos];
        else
            allSeeds = [allSeeds,cell2mat(tplan(i,2))'];
        end
    end
    
    % for correct data to be loaded, but patient id behind seed positions
    allSeeds = [allSeeds, patient];
    
    % evaluate gradients
    [fval,gradval] = matRad_calcGradients(allSeeds);
    doses = extractdata(fval);
    grads_x = extractdata(gradval(1));
    grads_y = extractdata(gradval(2));
    grads_z = extractdata(gradval(3));
    
    save(['..\BRACHYTHERAPY_data\evaluation\P_D90\fval_',num2str(patientId), '.mat'], "fval");
    save(['..\BRACHYTHERAPY_data\evaluation\P_D90\gradval_',num2str(patientId), '.mat'], "gradval");      
end