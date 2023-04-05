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
    [ct, cst, pln, stf] = matRad_automaticDifferentiationInit();

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
    refVol = 90;
    qi_all = matRad_calcQualityIndicators(cst,pln,doseCube,refGy,refVol);

    % Define wanted QI here
    qi = 100*qi_all(1,1).V_160Gy;
    dqidx = dlgradient(qi,x);

end

