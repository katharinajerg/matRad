function [ct, cst, pln, stf] = matRad_automaticDifferentiationInit(id)
%MATRAD_AUTOMATICDIFFERENTIATIONINIT is a matRad function to initialize all
% data for the differentiation of a qualitiy measure.
% 
% output
%   ct:         ct cube
%   cst:        matRad cst struct
%               (positions and constraints of patient structures)
%   pln:        matRad plan meta information struct
%   stf:        struct containing geometric information
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


%% I Patient Data Import
% Let's begin with a clear Matlab environment. Then, import the TG119
% phantom into your workspace. The phantom is comprised of a 'ct' and 'cst'
% structure defining the CT images and the structure set. Make sure the 
% matRad root directory with all its subdirectories is added to the Matlab 
% search path.
matRad_rc;
patient = id;
patientId = 1000+patient;
path = ['..\BRACHYTHERAPY_data\',num2str(patient),'\IntraOp\IntraOp\'];
pathStructureSet = [path, 'SS001.dcm']; 
pathImg = [path, 'MR001.dcm'];
pathPln = [path,'PL001.dcm'];

%% Import data
[cst, ct] = matRad_importDicomUSStructureSet(pathStructureSet,pathImg);
infoPl = dicominfo(pathPln);
targetDose = infoPl.DoseReferenceSequence.Item_1.TargetPrescriptionDose;


%% I - set dose objectives for brachytherapy
% Set the prostate as the Target amd remaining structures as OARs
% Prostate bed objective
cst{1,3} = 'TARGET';
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(400,targetDose));
cst{1,5}.Priority = 3;

% Rectum Objective
cst{3,3}    =  'OAR';
cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(30,targetDose));
cst{3,5}.Priority = 1;

% Urethra Objective
cst{2,3}    =  'OAR';
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(1,targetDose));
cst{2,5}.Priority = 2;

ct.info = 'automaticDiffInit';

save(['cst_',num2str(patientId), '.mat'], "cst");
save(['ct_',num2str(patientId), '.mat'], "ct");



%% II.1 Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% matlab structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon 
% (external beam) or brachy as an invasive tratment option.
% In this case we want to use brachytherapy. Then, we need to define a 
% treatment machine to correctly load the corresponding base data.
% matRad includes example base data for HDR and LDR brachytherapy.
% Here we will use HDR. By this means matRad will look for 'brachy_HDR.mat'
% in our root directory and will use the data provided in there for 
% dose calculation.

pln.radiationMode   = 'brachy'; 
pln.machine         = 'LDR';    % 'LDR' or 'HDR' for brachy


%% II.1 - needle and template geometry
% Now we have to set some parameters for the template and the needles. 
% Let's start with the needles: Seed distance is the distance between
% two neighbouring seeds or holding points on one needle or catheter. The
% seeds No denotes how many seeds/holding points there are per needle.

pln.propStf.needle.seedDistance      = 10; % [mm] 
pln.propStf.needle.seedsNo           = 6; 


%% II.1 - template position
% The implantation is normally done through a 13 x 13 template from the 
% patients inferior, which is the negative z axis here.
% The direction of the needles is defined by template normal.
% Neighbour distances are called by bixelWidth, because this field is also
% used for external beam therapy.
% The needles will be positioned right under the target volume pointing up.

pln.propStf.template.normal      = [0,0,1];
pln.propStf.bixelWidth   = 5; % [mm] template grid distance
%pln.propStf.templateRoot = matRad_getTemplateRoot(ct,cst); % mass center of
% target in x and y and bottom in z

% Here, we define active needles as 1 and inactive needles
% as 0. This is the x-y plane and needles point in z direction. 
% A checkerboard pattern is frequantly used. The whole geometry will become
% clearer when it is displayed in 3D view in the next section.

pln.propStf.template.activeNeedles = [0 0 0 1 0 1 0 1 0 1 0 0 0;... % 7.0
                                      0 0 1 0 1 0 0 0 1 0 1 0 0;... % 6.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 6.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 5.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 5.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 4.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 3.0
                                      1 0 1 0 1 0 1 0 1 0 1 0 1;... % 2.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 2.0
                                      1 0 1 0 1 0 0 0 0 0 1 0 1;... % 1.5
                                      0 0 0 0 0 0 0 0 0 0 0 0 0];   % 1.0
                                     %A a B b C c D d E e F f G

pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0); %  target center



%% II.1 - dose calculation options
% for dose calculation we use eather the 2D or the 1D formalism proposed by
% TG 43. Also, set resolution of dose calculation and optimization.
% If your system gets stuck with the resolution, you can lower it to 10 or
% 20, just to get an initial result. Otherwise, reduce the number of
% needles.
% Calculation time will be reduced by one tenth when we define a dose
% cutoff distance.


pln.propDoseCalc.TG43approximation = '1D'; %'1D' or '2D' 

% set dose cube resolution to same as ct resolution in order to avoid
% interpolation. Without optimization the size of the dose cube can as well
% be large.
pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z; % [mm]

pln.propDoseCalc.DistanceCutoff    = 130; %[mm] sets the maximum distance
                                            %to which dose is calculated. 

% the standard interior point optimizer IPOPT can be used
pln.propOpt.optimizer       = 'SA';

%% II.1 - book keeping
% Some field names have to be kept although they don't have a direct
% relevance for brachy therapy.
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runDAO          = false;  
pln.propOpt.runSequencing   = false; 
pln.propStf.gantryAngles    = []; 
pln.propStf.couchAngles     = []; 
pln.propStf.numOfBeams      = 0;
pln.propStf.importSeedPos   = 3;    % This is where to set the seeds!!! For path check matRad_generateBrachyStf()
pln.numOfFractions          = 1; 

%% II.1 - view plan
% Et voila! Our treatment plan structure is ready. Lets have a look:
disp(pln);


%% II.2 Steering Seed Positions From STF
% The steering file struct contains all needls/catheter geometry with the
% target volume, number of needles, seeds and the positions of all needles
% The one in the end enables visualization.
stf = matRad_generateStf(ct,cst,pln,1);
save(['stf_',num2str(patientId), '.mat'], "stf");


%% II.2 - view stf
% The 3D view is interesting, but we also want to know how the stf struct
% looks like.
disp(stf)

end

