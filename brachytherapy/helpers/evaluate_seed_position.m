% 14.03.2023, Katharina Jerg
% This script is used to deterimine the impact of each seed on different
% dose measures. 

clear all
close all

% 1. 
% in matRad_automaticDifferentiationInit.m:
% - l.78 set pln.machine = 'LDR';
% - l.137 set pln.propDoseCalc.TG43approximation = '1D';
% - l.158 set pln.propStf.importSeedPos = 3; This ensures that only the seeds of the input are taken for
%   differentiation. 

% 2. in matRad_automaticDifferentiation.m
% - l.15ff: define wanted quality indicators

% 3. import seed positions corresponding to the patient number
patients_160_dose_objectives = [1,4,10,11,41,42,51,57,59,62,64,70,72,87,90,96,99,101,109,117,120]; % patients with fullfilled 160Gy constraints
patients_110 = [10,11,16,35,64,66,72,78,90,95]; % patients with target dose 110Gy
patients_wrong_dose_objectives = [2,3,5,6,7,8,9,12,13,14,15,17,18,19,20,21,22,23,24,25]; % patients, which do not fullfill all 160Gy constraints
for patient = 70
    patientId = 1000+patient;
    load(['~/thindrives/ProstateData/',num2str(patient),'/IntraOp/IntraOp/tplan_orig.mat']);
    %load(['~/thindrives/ProstateData/Pat',num2str(patient),'/tplan_orig.mat']);
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
    
    save(['fval_',num2str(patientId), '.mat'], "fval");
    save(['gradval_',num2str(patientId), '.mat'], "gradval");
    
    % remove patient id from seeds again
    allSeeds = allSeeds(1:end-1);
    
    %%  visualization
    load(['gradval_',num2str(patientId) , '.mat']);
    load(['stf_',num2str(patientId) , '.mat']);
    gradPerSeed = reshape(extractdata(gradval), [3,numSeeds]);
    magnitudePerSeed = vecnorm(gradPerSeed);
    
    % seeds
    figure;
    set(gcf,'Position',[100 100 500 400]);
    hold on 
    S = 50;
    reshaped = reshape(allSeeds,3,numSeeds);
    for i = 1:size(tplan,1)
         p = scatter3(reshaped(1,i),reshaped(2,i),reshaped(3,i),S, magnitudePerSeed(i), 'filled');
    %      p = scatter3(tplan{i,2}(1), tplan{i,2}(2), tplan{i,2}(3), S, magnitudePerSeed(i), 'filled');
    end
    set(gca,'ColorScale','log')
    colorbar
    
    % prostate surface
    load(['ct_',num2str(patientId) , '.mat']);
    load(['cst_',num2str(patientId) , '.mat']);
    b = extract_tissue_boundaries(ct, cst);

    index_prostate = find(strcmp(b, 'Prostate'));
    index_urethra = find(strcmp(b, 'Urethra'));
    index_rectum = find(strcmp(b, 'Rectum'));
    
    s_prostate = b{index_prostate,2};
    P_prostate = b{index_prostate,3};
    s_urethra = b{index_urethra,2};
    P_urethra = b{index_urethra,3};
    s_rectum = b{index_rectum,2};
    P_rectum = b{index_rectum,3};

    trisurf(s_prostate,P_prostate(:,1),P_prostate(:,2),P_prostate(:,3),'FaceColor','red','FaceAlpha',0.1,'LineStyle','none')
    trisurf(s_urethra,P_urethra(:,1),P_urethra(:,2),P_urethra(:,3),'FaceColor','green','FaceAlpha',0.1,'LineStyle','none')
    trisurf(s_rectum,P_rectum(:,1),P_rectum(:,2),P_rectum(:,3),'FaceColor','blue','FaceAlpha',0.1,'LineStyle','none')
    xlabel('x')
    ylabel('y')
    zlabel('z')

    % brown:[0.8500 0.50 0.0980]
    view(40,30)
    
    % prostate contours
%     structureSetPath = ['~/thindrives/ProstateData/',num2str(patient),'/IntraOp/IntraOp/SS001.dcm'];
%     info = dicominfo(structureSetPath);
%     contour = dicomContours(info);
%     plotContour(contour)
    set(gcf,'PaperSize',[15 10]); %set the paper size to what you want  
    print(gcf,['~/thindrives/ProstateData/',num2str(patient),'/IntraOp/IntraOp/R_2cc.pdf'],'-dpdf')
    %saveas(gcf,['~/thindrives/ProstateData/',num2str(patient),'/IntraOp/IntraOp/V160.pdf'])
       
    %% save gif
    angleRange = 0:6:359;
    for i = 1:length(angleRange)
        % Rotate the camera view
        view(angleRange(i), 30);
        axis vis3d
        set(gcf,'Position',[100 100 500 400]);
        exportgraphics(figure(1),['~/thindrives/ProstateData/',num2str(patient),'/IntraOp/IntraOp/animation.gif'],'Append',true);
    end
    
       
    %%
    %%%%%%%%%% This is the code to move seeds of entire needles %%%%%%%%%%%%%%%%%% 
    
    % split seeds in needles
    % load(['gradval_',num2str(patientId),'.mat']);
    % gradPerSeed = reshape(extractdata(gradval), [3,size(gradval,2)/3]);
    % magnitudePerSeed = vecnorm(gradPerSeed);
    % numNeedles = 0;
    % seedIDsPerNeedle = cell(1);
    % IDs = 1;
    % for i = 2:numSeeds
    %     if(all(tplan{i,3} == [0;0;1]))
    %         if (tplan{i,2}(1) == tplan{i-1,2}(1))
    %             IDs = [IDs,i];
    %         else
    %             numNeedles = numNeedles +1;
    %             seedIDsPerNeedle{numNeedles} = IDs;
    %             IDs = i; 
    %         end
    %     else 
    %         warning('needle is not in z direction');
    %     end
    % end
    % numNeedles = numNeedles +1;
    % seedIDsPerNeedle{numNeedles} = IDs;
    % 
    % % find mean gradient magnitude per needle
    % magnitudePerNeedle = zeros(1,numNeedles);
    % for i = 1:numNeedles
    %     IDs = seedIDsPerNeedle{i};
    %     magnitudePerNeedle(i) = mean(magnitudePerSeed(IDs));
    % end
    % 
    % % move needles dependent on gradient
    % [maxGrad, maxNeedleId] = max(magnitudePerNeedle);
    % maxSeedIds = seedIDsPerNeedle{maxNeedleId};


end