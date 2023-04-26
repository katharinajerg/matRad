function [qi, dqidx] = matRad_automaticDifferentiation(x)
    %MATRAD_AUTOMATICDIFFERENTIATION is a matRad function to differentiate 
    % a qualitiy indicator depending on given seed positions. 
    % 
    % input
    % x:        Vector of seed coordinates [x1,y1,z1, ..., xn, yn, zn] for 
    %           n seeds
    % 
    % output
    % qi:       function value of the desired quality measure
    % dqidx:    derivative of the quality measure regarding each input
    %           coordinate
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2021 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Init structs
    id = extractdata(x(end));
    x = x(1:end-1);

    [ct, cst, pln, stf] = matRad_automaticDifferentiationInit(id);


    % check for seeds exactly on voxel boundary to avoid NaN gradients
    seeds = reshape(extractdata(x), [3,numel(x)/3]);
    for s = 1:size(seeds,2)
        if( ismember(seeds(1,s), ct.x) && ismember(seeds(2,s), ct.y) && ismember(seeds(3,s), ct.z))
            x((s-1)*3+1) = x((s-1)*3+1)+0.0001;
        end
    end
   
    % calculate dose influence matrix
    dij = matRad_calcBrachyDose(ct,stf,pln,cst,x);
 
    % calculate dose cube from dose influence matrix
    w = ones(numel(x)/3,1);
    resultGUI = matRad_calcCubes(w,dij);
    doseCube = resultGUI.physicalDose;

    % plot dose and ct slice for visualization
    doseSlice = extractdata(doseCube(:,:,floor(size(doseCube,3)/2)));
    ctCube = ct.cube{1,1};
    ctSlice = ctCube(:,:,floor(size(ctCube,3)/2));
    ctSliceFiltered = edge(ctSlice);
    doseSlice = doseSlice + max(doseSlice(:))*ctSliceFiltered;
    figure
    imagesc(doseSlice)
    
    % calculate QI
    refGy = 160;
    refVol = 30;
    qi_all = matRad_calcQualityIndicators(cst,pln,doseCube,refGy,refVol);


    qi = qi_all(3).D_2CC;
    dqidx = dlgradient(qi,x);

    % dose contraints
    targetDose = 160;
    [dvh,qi_2] = matRad_indicatorWrapper(cst,pln,resultGUI, [targetDose,1.5*targetDose,2*targetDose], ...
    [10, 30, 90, 95]);
    fprintf(extractdata(1+~(qi_2(1).D_90 > (0.9*targetDose))), 'P_D90 is %5.2fGy and should be greater than %5.2fGy\n', qi_2(1).D_90, 0.9*targetDose);
    fprintf(extractdata(1+~(eval(['qi_2(1).V_',num2str(targetDose),'Gy']) > 0.95)), 'P_V100 is %5.2f and should be greater than 0.95\n', eval(['qi_2(1).V_',num2str(targetDose),'Gy']));
    fprintf(extractdata(1+~(eval(['qi_2(1).V_',num2str(1.5*targetDose),'Gy']) > 0.45)+~(eval(['qi_2(1).V_',num2str(1.5*targetDose),'Gy']) < 0.65)), 'P_V150 is %5.2f and should be between 0.45 and 0.65\n', eval(['qi_2(1).V_',num2str(1.5*targetDose),'Gy']));
    fprintf(extractdata(1+~(qi_2(2).D_10 < 240)), 'U_D10 is %5.2fGy and should be less than 240Gy\n', qi_2(2).D_10);
    fprintf(extractdata(1+~(qi_2(2).D_30 < 208)), 'U_D30 is %5.2fGy and should be less than 208Gy\n', qi_2(2).D_30);
    fprintf(extractdata(1+~(qi_2(3).D_01CC < 200)), 'R_D01CC is %5.2fGy and should be less than 200Gy\n', qi_2(3).D_01CC);
    if(isfield(qi_2, 'D_2CC'))
        fprintf(extractdata(1+~(qi_2(3).D_2CC < 145)), 'R_D2CC is %5.2fGy and should be less than 145Gy\n', qi_2(3).D_2CC);
    end
end

