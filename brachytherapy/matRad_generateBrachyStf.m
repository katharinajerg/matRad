function stf = matRad_generateBrachyStf(ct,cst,pln, visMode)
% matRad_generateBrachyStf generates matRad steering information generation for brachy
%   will be called within matRad_generateStf if radiation mode is 'brachy'
%
% call
%   stf = matRad_generateBrachyStf(ct,cst,pln,visMode) 
%   
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%   pln:        matRad plan meta information struct
%   visMode:    toggle on/off different visualizations by setting this 
%               value to 1,2,3 (optional)
%
% output
%   stf:        matRad steering information struct
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

%% config
matRad_cfg = MatRad_Config.instance();
addpath(fullfile( matRad_cfg.matRadRoot));
matRad_cfg.dispInfo('matRad: Generating stf struct... ');

if ~isfield(pln,'propStf')
matRad_cfg.dispError('no applicator information in pln struct');
end



%% generate image coordinates

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
    
% add margin
addmarginBool = matRad_cfg.propStf.defaultAddMargin;
if isfield(pln,'propStf') && isfield(pln.propStf,'addMargin')
   addmarginBool = pln.propStf.addMargin; 
end

if addmarginBool
    voiTarget = matRad_addMargin(voiTarget,cst,ct.resolution,ct.resolution,true);
    V   = find(voiTarget>0);
end

% throw error message if no target is found
if isempty(V)
    matRad_cfg.dispError('Could not find target.');
end

% Convert linear indices to 3D voxel coordinates
[coordsY_vox, coordsX_vox, coordsZ_vox] = ind2sub(ct.cubeDim,V);

% % calculate rED or rSP from HU
% ct = matRad_calcWaterEqD(ct, pln);

% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));

% ignore densities outside of contours
eraseCtDensMask = ones(prod(ct.cubeDim),1);
eraseCtDensMask(V) = 0;
for i = 1:ct.numOfCtScen
    ct.cube{i}(eraseCtDensMask == 1) = 0;
end

%% save VOI coordinates
%translate to geometric coordinates and save in stf
stf.targetVolume.Xvox = ct.x(coordsX_vox);
stf.targetVolume.Yvox = ct.y(coordsY_vox);
stf.targetVolume.Zvox = ct.z(coordsZ_vox);


 %% meta info from pln
 stf.radiationMode = pln.radiationMode;
   
        
%% generate seed positions
% seed positions can be generated from needles, template and orientation
% needles are assumed to go through the template vertically

% needle position
% when import does not exist or is false create seed positons from template
if (~isfield(pln.propStf, 'importSeedPos')| ~pln.propStf.importSeedPos)

    %% generate 2D template points
    % the template origin is set at its center. In the image coordinate system,
    % the center will be positioned at the bottom of the volume of
    % interest.
    [row,col] = find(pln.propStf.template.activeNeedles);
    templX = col*pln.propStf.bixelWidth + pln.propStf.templateRoot(1) - (13+1)/2*pln.propStf.bixelWidth;
    templY = row*pln.propStf.bixelWidth + pln.propStf.templateRoot(2) - (13+1)/2*pln.propStf.bixelWidth;
    templZ = ones(size(col))                 + pln.propStf.templateRoot(3);
    
    stf.template = [templX';templY';templZ'];

    % more meta data
    stf.numOfSeedsPerNeedle = pln.propStf.needle.seedsNo;
    stf.numOfNeedles = nnz(pln.propStf.template.activeNeedles);
    stf.totalNumOfBixels = stf.numOfSeedsPerNeedle*stf.numOfNeedles; % means total number of seeds 

    d = pln.propStf.needle.seedDistance;
    seedsNo = pln.propStf.needle.seedsNo;
    needleDist(1,1,:) = d.*[0:seedsNo-1]'; % 1x1xN Array with seed positions on needle
    needleDir = needleDist.*[0;0;1];
    seedPos_coord_need_seed = needleDir + stf.template;
    seedPos_need_seed_coord = shiftdim(seedPos_coord_need_seed,1);
    % the output array has the dimensions (needleNo,seedNo,coordinates)

    X = seedPos_need_seed_coord(:,:,1);
    Y = seedPos_need_seed_coord(:,:,2);
    Z = seedPos_need_seed_coord(:,:,3);
    
    stf.seedPoints.x = reshape(X,1,[]);
    stf.seedPoints.y = reshape(Y,1,[]);
    stf.seedPoints.z = reshape(Z,1,[]);

    seedpointSize = size(stf.seedPoints.x);
    stf.seedPoints.x_orientation = zeros(seedpointSize);
    stf.seedPoints.y_orientation = zeros(seedpointSize);
    stf.seedPoints.z_orientation = ones(seedpointSize);

elseif (pln.propStf.importSeedPos == 1) 
    % load dwell points from tplan file
    % tplan = calcDwellPoints()
    
    input = load('./brachytherapy/data/tplan_full.mat', 'full_tplan');
    full_tplan = input.full_tplan;
    
    x = zeros(1,size(full_tplan,1));
    y = zeros(1,size(full_tplan,1));
    z = zeros(1,size(full_tplan,1));

    x_orientation = zeros(1,size(full_tplan,1));
    y_orientation = zeros(1,size(full_tplan,1));
    z_orientation = zeros(1,size(full_tplan,1));

    for i = 1:size(full_tplan,1)
        x(i) = full_tplan{i,2}(1);
        y(i) = full_tplan{i,2}(2);
        z(i) = full_tplan{i,2}(3);

        x_orientation(i) = full_tplan{i,3}(1);
        y_orientation(i) = full_tplan{i,3}(1);
        z_orientation(i) = full_tplan{i,3}(1);
    end

    stf.seedPoints.x = x;
    stf.seedPoints.y = y;
    stf.seedPoints.z = z;

    stf.seedPoints.x_orientation = x_orientation;
    stf.seedPoints.y_orientation = y_orientation;
    stf.seedPoints.z_orientation = z_orientation;

    % more meta data
    stf.numOfSeedsPerNeedle = [];
    stf.numOfNeedles = [];
    stf.totalNumOfBixels = numel(stf.seedPoints.x); % means total number of seeds 
    stf.template = [];

elseif (pln.propStf.importSeedPos == 2)
    %load support points, which represent the needles' paths from file

    load('./brachytherapy/data/suppPoints1.mat', 'supportPoints');
    numNeedles = size(supportPoints,1);
    needleController = MatRad_BrachyGeometryController(supportPoints, numNeedles, ...
        pln.propStf.needle.seedDistance, pln.propStf.needle.seedsNo);
    needleController = needleController.calcNeedles();
    
    x = [];
    y = [];
    z = [];
    x_orientation = [];
    y_orientation = [];
    z_orientation = [];

    for n = 1:needleController.numberOfNeedles
        for d = 1:size(needleController.needleSet{n}.dwellPointList,2)
            x  = [x, needleController.needleSet{n}.dwellPointList{d}.positionX];
            y  = [y, needleController.needleSet{n}.dwellPointList{d}.positionY];
            z  = [z, needleController.needleSet{n}.dwellPointList{d}.positionZ];

            x_orientation  = [x_orientation, needleController.needleSet{n}.dwellPointList{d}.orientationX];
            y_orientation  = [y_orientation, needleController.needleSet{n}.dwellPointList{d}.orientationY];
            z_orientation  = [z_orientation, needleController.needleSet{n}.dwellPointList{d}.orientationZ];
        end
    end

    stf.seedPoints.x = x;
    stf.seedPoints.y = y;
    stf.seedPoints.z = z;

    stf.seedPoints.x_orientation = x_orientation;
    stf.seedPoints.y_orientation = y_orientation;
    stf.seedPoints.z_orientation = z_orientation;

    % more meta data
    stf.numOfSeedsPerNeedle = pln.propStf.needle.seedsNo;
    stf.numOfNeedles = needleController.numberOfNeedles;
    stf.totalNumOfBixels = numel(stf.seedPoints.x); % means total number of seeds 
    stf.template = [];

elseif (pln.propStf.importSeedPos == 3)
    %do not define any seed positions. Will be defined later. Used for
    %automatic differentiation of seed position.

    stf.seedPoints.x = [];
    stf.seedPoints.y = [];
    stf.seedPoints.z = [];

    stf.seedPoints.x_orientation = [];
    stf.seedPoints.y_orientation = [];
    stf.seedPoints.z_orientation = [];

    % more meta data
    stf.numOfSeedsPerNeedle = [];
    stf.numOfNeedles = [];
    stf.totalNumOfBixels = numel(stf.seedPoints.x); % means total number of seeds 
    stf.template = [];

else
    error('No known action where to get the dwell points from. pln.propStf.importSeedPos must be either 0, 1, or 2.')
end

matRad_cfg.dispInfo('...100% ');

%%visualize results of visMode is nonzero
% plot 3D seed positions
if visMode > 0
    clf
    SeedPoints = plot3(stf.seedPoints.x,stf.seedPoints.y,stf.seedPoints.z,'.','DisplayName', 'seed points','Color','black','markersize',5);
    title( '3D Visualization of seed points')
    xlabel('X (left) [mm]')
    ylabel('Y (posterior) [mm]')
    zlabel('Z (superior) [mm]')
    hold on
    

    % plot 3d VOI points
    TargX = stf.targetVolume.Xvox;
    TargY = stf.targetVolume.Yvox;
    TargZ = stf.targetVolume.Zvox;
    %Prostate = plot3(TargX,TargY,TargZ,'.', 'Color','b','DisplayName', 'prostate');
    
    P = [TargX',TargY',TargZ'];
    k = boundary(P,1);
    trisurf(k,P(:,1),P(:,2),P(:,3),'FaceColor','red','FaceAlpha',0.1,'LineStyle','none')
    hold off;
end

end

