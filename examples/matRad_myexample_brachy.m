clear all 

matRad_rc;
pathStructureSet = "~/Daten/SS001.dcm"; 
pathImg = "~/Daten/MR001.dcm";
[cst, ct] = matRad_importDicomUSStructureSet(pathStructureSet,pathImg);

%% I - set dose objectives for brachytherapy

% Set the prostate as the Target amd remaining structures as OARs

% Prostate bed objective
cst{1,3} = 'TARGET';
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(1,140));
cst{1,5}.Priority = 3;

% Rectum Objective
cst{3,3}    =  'OAR';
cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,120));
cst{3,5}.Priority = 1;

% Urethra Objective
cst{2,3}    =  'OAR';
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,120));
cst{2,5}.Priority = 2;


%% II - Set seed geometry and planning parameters and calculate the dose influence matrix

% II.1 Treatment Plan
pln.radiationMode   = 'brachy'; 
pln.machine         = 'LDR';    % 'LDR' or 'HDR' for brachy

% II.2 - needle and template geometry
pln.propStf.importSeedPos           = 1; % 1 for true (seed positions are imported from tplan.mat file), 0 for false 
pln.propStf.needle.seedDistance     = 5; % [mm] 
pln.propStf.needle.seedsNo          = 5; 

% II.3 - template position
pln.propStf.template.normal = [0,0,1];
pln.propStf.bixelWidth      = 5; % [mm] template grid distance
pln.propStf.templateRoot    = matRad_getTemplateRoot(ct,cst); 
% mass center of target in x and y and bottom in z

% Here, we define active needles as 1 and inactive needles
% as 0. This is the x-y plane and needles point in z direction. 
pln.propStf.template.activeNeedles = [0 0 0 0 0 0 0 0 0 0 0 0 0;... % 7.0
                                      0 0 0 0 0 0 0 0 0 0 0 0 0;... % 6.5
                                      0 0 0 0 0 0 0 0 0 0 0 0 0;... % 6.0
                                      0 0 0 0 0 0 0 0 0 0 0 0 0;... % 5.5
                                      0 0 0 0 1 1 1 1 1 0 0 0 0;... % 5.0
                                      0 0 1 1 1 1 1 1 1 1 1 0 0;... % 4.5
                                      0 0 1 1 1 1 1 1 1 1 1 0 0;... % 4.0
                                      0 0 1 1 1 1 1 1 1 1 1 0 0;... % 3.5
                                      0 0 0 0 1 1 1 1 1 0 0 0 0;... % 3.0
                                      0 0 0 0 0 0 0 0 0 0 0 0 0;... % 2.5
                                      0 0 0 0 0 0 0 0 0 0 0 0 0;... % 2.0
                                      0 0 0 0 0 0 0 0 0 0 0 0 0;... % 1.5
                                      0 0 0 0 0 0 0 0 0 0 0 0 0];   % 1.0
                                     %A a B b C c D d E e F f G

pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0); %  target center

% II.4 - dose calculation options
pln.propDoseCalc.TG43approximation = '1D'; %'1D' for LDR or '2D' for HDR 

pln.propDoseCalc.doseGrid.resolution.x = 1; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 1; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 2.5; % [mm]

pln.propDoseCalc.DistanceCutoff    = 130; % [mm] sets the maximum distance
                                          % to which dose is calculated. 

% the optimizer SA is used
pln.propOpt.optimizer       = 'SA';

% II.5 - book keeping
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runDAO          = false;  
pln.propOpt.runSequencing   = false; 
pln.propStf.gantryAngles    = []; 
pln.propStf.couchAngles     = []; 
pln.propStf.numOfBeams      = 0;
pln.numOfFractions          = 1; 

% II.6 - view plan
disp(pln);

% II.7 Steering Seed Positions From STF
figure
stf = matRad_generateStf(ct,cst,pln,1);

% II.8 - view stf
% The 3D view is interesting, but we also want to know how the stf struct
% looks like.
disp(stf)

% II.9 - Dose Calculation
dij = matRad_calcBrachyDose(ct,stf,pln,cst);
save("dij.mat", "dij")


%% III - Inverse Optimization for brachy therapy

resultGUI = matRad_fluenceOptimization(dij,cst,pln);
save("result.mat", "resultGUI")
save("ct.mat", "ct")
matRadGUI;


%% IV - Calculate and plot the results

% IV.1 plot the transversal iso-center dose slice
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
figure
imagesc(resultGUI.physicalDose(:,:,slice)),colorbar, colormap(jet);

% IV.2 Obtain dose statistics
% Two more columns will be added to the cst structure depicting the DVH and
% standard dose statistics such as D95,D98, mean dose, max dose etc.
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI, [140,210,280], [5, 90, 95, 98]);


