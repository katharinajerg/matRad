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

    % calculate dose influence matrix
    dij = matRad_calcBrachyDose(ct,stf,pln,cst,x);
 
    % calculate dose cube from dose influence matrix
    w = ones(numel(x)/3,1);
    resultGUI = matRad_calcCubes(w,dij);
    doseCube = resultGUI.physicalDose;

    % calculate QI
    refGy = 160;
    refVol = 90;
    qi_all = matRad_calcQualityIndicators(cst,pln,doseCube,refGy,refVol);

    % !!! Define wanted QI here!!!
%     qiD = qi_all(1,1).D_90;
    qi = qi_all(1,1).V_160Gy;
    %dqiddx = dlgradient(qiD,x);
    dqidx = dlgradient(qiV,x);

end

