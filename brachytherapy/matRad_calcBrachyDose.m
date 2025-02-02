function dij = matRad_calcBrachyDose(ct,stf,pln,cst,x,id)
% matRad_calcBrachyDose calculates dose influence matrix according to the
% AAPM update Rivard et al. 2004
%
% call
%   dij = matRad_calcBrachyDose(ct,stf,pln,cst)
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%               (positions and constraints of patient structures)
%   pln:        matRad plan meta information struct
%   stf:        struct containing geometric information
%   x:          (optional) 3n vector with seed positions (x_1, y_1, z_1, ..., y_n, z_n)
%   id:         (optional) patient id for loading machine data
%
% output
%   dij:        stuct containing dose influence information
%
% References: 
%   [1] https://doi.org/10.1118/1.1646040 - TG43 Update
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

if ~exist('id','var') || isempty(id) 
    id = 1;
end

%%Configure
matRad_cfg =  MatRad_Config.instance();
matRad_calcBrachyDoseInit;

% initialize waitbar (always indented to seperate from important code)
figureWait = waitbar...
    (0,'calculating dose inlfluence matrix for brachytherapy...');
matRad_cfg.dispInfo('Starting  brachytherapy dose calculation...\n');
startTime = tic;

%% get dose points and seedpoints
% "dosePoints" and "seedPoints" are both structs with fields x,y,z:
% each contains a 1D row vector of position components [mm]

if ~exist('x','var') || isempty(x) 
    seedPoints.x = single(stf.seedPoints.x);
    seedPoints.y = single(stf.seedPoints.y);
    seedPoints.z = single(stf.seedPoints.z);
    
    seedPoints.x_orientation = single(stf.seedPoints.x_orientation);
    seedPoints.y_orientation = single(stf.seedPoints.y_orientation);
    seedPoints.z_orientation = single(stf.seedPoints.z_orientation);
else
    numSeedCoordinates = numel(x);
    assert(mod(numSeedCoordinates, 3)==0, 'The size of the vector defining the seed positions can not be devided by 3. Check seed position vector.\n')
    numSeeds = numSeedCoordinates/3;
    seeds = reshape(x,[3,numSeeds]);

    seedPoints.x = single(seeds(1,:));
    seedPoints.y = single(seeds(2,:));
    seedPoints.z = single(seeds(3,:));

    seedPoints.x_orientation = zeros(1,numSeeds);
    seedPoints.y_orientation = zeros(1,numSeeds);
    seedPoints.z_orientation = ones(1,numSeeds);
end



[XGrid,YGrid,ZGrid] = meshgrid(dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
dosePoints.x = single(reshape(XGrid,1,[]));
dosePoints.y = single(reshape(YGrid,1,[]));
dosePoints.z = single(reshape(ZGrid,1,[]));

matRad_cfg.dispInfo('\t computing distance transform... ');


%% get seed dosepoint distance matrix
% [seedPoint x dosePoint] matrix with relative distance as entries
% detailed documentation in function
DistanceMatrix = matRad_getDistanceMatrix(seedPoints,dosePoints);
OrientationMatrix =  matRad_getOrientationMatrix(seedPoints,dosePoints);

% ignore all distances > Cutoff for the following calculations to save time
Ignore = DistanceMatrix.dist > pln.propDoseCalc.DistanceCutoff;
calcDistanceMatrix.x = DistanceMatrix.x(~Ignore);
calcDistanceMatrix.y = DistanceMatrix.y(~Ignore);
calcDistanceMatrix.z = DistanceMatrix.z(~Ignore);
calcDistanceMatrix.dist = DistanceMatrix.dist(~Ignore);

calcOrientationMatrix.x = OrientationMatrix.x_orientation(~Ignore);
calcOrientationMatrix.y = OrientationMatrix.y_orientation(~Ignore);
calcOrientationMatrix.z = OrientationMatrix.z_orientation(~Ignore);
% remove singularities
calcDistanceMatrix.dist(calcDistanceMatrix.dist < machine.data.ActiveSourceLength) = machine.data.ActiveSourceLength;

% now all fields of calcDistanceMatrix are n x 1 arrays!

% update waitbar
waitbar(0.125);
matRad_cfg.dispInfo('done in %f s!\n',toc(startTime));

%% seed dosepoint angle matrix
% [seedPoint x dosePoint] matrix with relative theta angle as entries
% detailed documentation in function
% only call for 2D formalism
if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'TG43approximation')
    pln.propDoseCalc.TG43approximation = '2D';
end

if strcmp(pln.propDoseCalc.TG43approximation,'2D')
    matRad_cfg.dispInfo('\t computing angle for TG43-2D... ');
    tmpTimer = tic;
    [ThetaMatrix,~] = matRad_getThetaMatrix(calcOrientationMatrix,calcDistanceMatrix);
    matRad_cfg.dispInfo('done in %f s!\n',toc(tmpTimer));
end

% update waitbar
waitbar(0.25);
%% Calculate Dose Rate matrix
% Calculation according to [1]

matRad_cfg.dispInfo('\t computing dose-rate for TG43-%s... ',pln.propDoseCalc.TG43approximation);
tmpTimer = tic;
if (exist('x','var') && isdlarray(x))
    DoseRate = dlarray(zeros(length(dosePoints.x),length(seedPoints.x)));
else
    DoseRate = zeros(length(dosePoints.x),length(seedPoints.x));
end

switch pln.propDoseCalc.TG43approximation
    case '1D'        
        DoseRate(~Ignore) = ...
        matRad_getDoseRate1D_poly(machine,calcDistanceMatrix.dist);
    case '2D'
        DoseRate(~Ignore) = ...
        matRad_getDoseRate2D_poly(machine,calcDistanceMatrix.dist,ThetaMatrix);
    otherwise
        matRad_cfg.dispError('TG43 Approximation ''%s'' not known!',pln.propDoseCalc.TG43approximation);
end
matRad_cfg.dispInfo('done in %f s!\n',toc(tmpTimer));

if (strcmp(pln.machine,'LDR'))
    % convert DoseRate in cGy per hour to Gy
    % integration over time from 0 -> inf yields
    dose = DoseRate*0.01*machine.data.SourceIsotopeHalfLife*24/0.693147; % /ln(2)
    dij.physicalDose = {dose}; % dose in Gy
else 
    dij.physicalDose = {DoseRate}; % dose rate in cGy per hours
end

% update waitbar, delete waitbar
waitbar(1);
matRad_cfg.dispInfo('Brachytherapy dose calculation finished in %f s!\n',toc(startTime));
delete(figureWait);


end



