% 14.03.2023, Katharina Jerg
% This script is used to deterimine the impact of each seed on different
% dose measures. 

clear all
close all

% 1. 
% in matRad_automaticDifferentiationInit.m:
% - l.11 set patient;
% - l.78 set pln.machine = 'LDR';
% - l.137 set pln.propDoseCalc.TG43approximation = '1D';
% - l.158 set pln.propStf.importSeedPos = 3; This ensures that only the seeds of the input are taken for
%   differentiation. 

% 2. in matRad_automaticDifferentiation.m
% - l.15ff: define wanted quality indicators

% 3. import seed positions corresponding to the patient number
% patient = 2;
%     load(['/home/kjerg/Results/2022_07_06 needle geometries/scan+plan/Pat', num2str(patient), '/tplan_orig.mat']);
%     numSeeds = size(tplan, 1);
%     allSeeds = [];
%     for i = 1:numSeeds
%         allSeeds = [allSeeds,cell2mat(tplan(i,2))'];
%     end

    tplan = cell(1,3);
    numSeeds = 1;
    allSeeds = [-2.01,0,0];
    tplan{1,2} = allSeeds';
    
    % evaluate gradients
    [fval,gradval] = matRad_calcGradients(allSeeds);
    doses = extractdata(fval);
    grads_x = extractdata(gradval(1));
    grads_y = extractdata(gradval(2));
    grads_z = extractdata(gradval(3));
    
%     save(['fval_',num2str(patient), '.mat'], "fval");
%     save(['gradval_',num2str(patient), '.mat'], "gradval");
 
%     %%   visualization
%     load(['gradval_',num2str(patient) , '.mat']);
%     load(['stf_',num2str(patient) , '.mat']);
%     gradPerSeed = reshape(extractdata(gradval), [3,numSeeds]);
%     magnitudePerSeed = vecnorm(gradPerSeed);
%     
%     % seeds
%     figure
%     hold on 
%     S = 50;
%     for i = 1:size(tplan,1)
%          p = scatter3(tplan{i,2}(1), tplan{i,2}(2), tplan{i,2}(3), S, magnitudePerSeed(i), 'filled');
%     end
%     set(gca,'ColorScale','log')
% 
%     % prostate surface
%     TargX = stf.targetVolume.Xvox;
%     TargY = stf.targetVolume.Yvox;
%     TargZ = stf.targetVolume.Zvox;
%     
%     P = [TargX',TargY',TargZ'];
%     k = boundary(P,1);
%     trisurf(k,P(:,1),P(:,2),P(:,3),'FaceColor','red','FaceAlpha',0.1,'LineStyle','none')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
% 
%     % prostate contours
% %     structureSetPath = ['~/Daten/Pat',num2str(patient),'/scan+plan/SS001.dcm'];
% %     info = dicominfo(structureSetPath);
% %     contour = dicomContours(info);
% %     plotContour(contour)
% 
%     view([40 10])
% 
%     
% %     %% save gif
% %     for n = 0:72
% %       view([n*5 10])
% %       exportgraphics(gcf,'Animated_2.gif','Append',true);
% %     end

