% This script is used to evaluate the results of the differentiation of the
% dose contraints regarding the seed positions.

%%
clear all
close all

% patient data sets
%fulfilling 160Gy dose constraints
% patientIds = [1001, 1004, 1010, 1011, 1041, 1042, 1051, 1057, 1059, 1062, 1064, 1070, 1072, 1087, 1090, 1096, 1099, 1101, 1109, 1117, 1120];
patientIds = 1004;% [1001:1009,1012:1015,1017:1025, 1041, 1042, 1051, 1057, 1059, 1062, 1064, 1070, 1087, 1096, 1099, 1101, 1109, 1117, 1120];
data = cell(size(patientIds,2),1);

% load data
for i = 1:size(patientIds,2)
    id = patientIds(i);
    load(['~/Results/2023_04_11_differentiation_on_patient_data/P_V150/gradval_',num2str(id),'.mat']);
    load(['~/Results/2023_04_11_differentiation_on_patient_data/P_V150/fval_',num2str(id),'.mat']);

    % V100 gradients 
    gradPerSeed = reshape(extractdata(gradval), [3,size(gradval,2)/3]);
    magnitudePerSeed = vecnorm(gradPerSeed);
    data{i}.id = id;
    data{i}.V100_gradients = gradPerSeed;
    data{i}.V100_gradMagnitudes = magnitudePerSeed;
    data{i}.V100 = extractdata(fval);
    data{i}.mean_gradMagnitudes = mean(magnitudePerSeed);
    data{i}.max_gradMagnitudes = max(magnitudePerSeed);

    % prostate information
    path = ['~/thindrives/ProstateData/',num2str(id-1000),'/IntraOp/IntraOp/'];
    pathStructureSet = [path, 'SS001.dcm'];
    pathImg = [path, 'MR001.dcm'];
    [cst, ct] = matRad_importDicomUSStructureSet(pathStructureSet,pathImg);
    cube = zeros(size(ct.cube{1}));
    cube(cst{1,4}{1}) = 1;
    data{i}.prostateVolume = sum(cube(:));

    %   extend in x, y, z
    line_x = sum(cube,3);
    line_x = sum(line_x,1); % x is dim = 2
    x_min =  ct.x(find(line_x, 1, 'first'));
    x_max =  ct.x(find(line_x, 1, 'last'));
    data{i}.xExtendProstate = [x_min, x_max];  

    line_y = sum(cube,3);
    line_y = sum(line_y,2); % y is dim = 1
    y_min =  ct.y(find(line_y, 1, 'first'));
    y_max =  ct.y(find(line_y, 1, 'last'));
    data{i}.yExtendProstate = [y_min, y_max];

    line_z = sum(cube,1);
    line_z = sum(line_z,2);
    z_min =  ct.z(find(line_z, 1, 'first'));
    z_max =  ct.z(find(line_z, 1, 'last'));
    data{i}.zExtendProstate = [z_min, z_max];

    %   center of mass
    [X, Y, Z] = meshgrid(ct.x, ct.y, ct.z);
    x_com = sum(sum(sum(X.*cube))) / sum(cube(:));
    y_com = sum(sum(sum(Y.*cube))) / sum(cube(:));
    z_com = sum(sum(sum(Z.*cube))) / sum(cube(:));
    data{i}.centerOfMass = [x_com, y_com, z_com];

    % seed positions
    load(['~/thindrives/ProstateData/',num2str(id-1000),'/IntraOp/IntraOp/tplan_orig.mat']);
    data{i}.tplan = tplan;

end


%% evaluate data
close all
% first overview over gradients
magnitudes = [];
g = [];
for i=1:size(patientIds,2)
    magnitudes = [magnitudes; data{i}.V100_gradMagnitudes'];
    g = [g; repmat({num2str(data{i}.id)}, size(data{i}.V100_gradMagnitudes,2),1)];

end

figure(1)
hold on
boxplot(magnitudes,g)
xlabel('patientID')
ylabel('V100 gradient magnitude')
title('first overview of V100 gradients')

% split data in "equally distributed" and "outlier seeds"
magnitudesEqual = [];
magnitudesOutlier = [];
gEqual = [];
gOutlier = [];
maxGrad = [];
meanGrad = [];
threshold = 3;
for i = 1:size(patientIds,2)
    if (max(data{i}.V100_gradMagnitudes) > threshold)
        magnitudesOutlier = [magnitudesOutlier; data{i}.V100_gradMagnitudes'];
        gOutlier = [gOutlier; repmat({num2str(data{i}.id)}, size(data{i}.V100_gradMagnitudes,2),1)];
        maxGrad = [maxGrad, max(data{i}.V100_gradMagnitudes)];
        meanGrad = [meanGrad, median(data{i}.V100_gradMagnitudes)];
    else
        magnitudesEqual = [magnitudesEqual; data{i}.V100_gradMagnitudes'];
        gEqual = [gEqual; repmat({num2str(data{i}.id)}, size(data{i}.V100_gradMagnitudes,2),1)];
        maxGrad = [maxGrad, max(data{i}.V100_gradMagnitudes)];
        meanGrad = [meanGrad, median(data{i}.V100_gradMagnitudes)];
    end
end

figure(2) 
boxplot(magnitudesEqual,gEqual)
xlabel('patientID')
ylabel('V100 gradient magnitude')
title('V100 gradients with max gradient < 0.6')

figure(3) 
boxplot(magnitudesOutlier,gOutlier)
xlabel('patientID')
ylabel('V100 gradient magnitude')
title('V100 gradients with max gradient > 0.6')

figure(4)
hold on
[sortedMaxs, ind] = sort(maxGrad);
sortedIds = patientIds(ind);
plot(1:size(patientIds,2), sort(maxGrad))
plot(1:size(patientIds,2),linspace(threshold,threshold + 0.0001,size(patientIds,2)))
hold off
xlabel('patientID')
ylabel('max(V100 gradient magnitude)')

% -->> The gradients of V100 per seed differ much between patients. The max
%      gradient per patient is between 0.047 and 8.26.

%% understand orientation
gradX = [];
gradY = [];
gradZ = [];
g = [];
for i = 1:size(patientIds,2)
    gradX = [gradX; data{i}.V100_gradients(1,:)'];
    gradY = [gradY; data{i}.V100_gradients(2,:)'];
    gradZ = [gradZ; data{i}.V100_gradients(3,:)'];
    g = [g; repmat({num2str(data{i}.id)}, size(data{i}.V100_gradMagnitudes,2),1)];
end
gradX = sort(abs(gradX));
gradY = sort(abs(gradY));
gradZ = sort(abs(gradZ));

figure(5)
plot(1:size(gradX,1), gradX, 1:size(gradY,1), gradY, 1:size(gradZ,1), gradZ);
xlabel('seed')
ylabel('V100 gradient')
title('V100 gradients by direction')
legend('x', 'y', 'z')

% -->> The gradients are mostly homogeneous for different orientations.

%% dependence on other variables
% V100
V100 = [];
for i = 1:size(patientIds,2)
    V100 = [V100,data{i}.V100];
end
figure(6)
hold on
yyaxis left
scatter(V100, maxGrad, 'blue', 'filled')
xlabel('V100 in %')
ylabel('max gradient magnitude')
yyaxis right
scatter(V100, meanGrad, 'red', 'filled')
ylabel('mean gradient magnitude')
hold off

% -->> There is no clear linear dependence from the gradient to the V100
%      value. There is a trend that for larger V100, the gradients of
%      single seeds is larger.


% Number of Seeds
numSeeds = [];
for i = 1:size(patientIds,2)
    numSeeds = [numSeeds,size(data{i}.V100_gradients,2)];
end
figure(7)
hold on
yyaxis left
scatter(numSeeds, maxGrad,'blue', 'filled')
xlabel('number of seeds per patient')
ylabel('max gradient magnitude')
yyaxis right 
scatter(numSeeds, meanGrad,'red', 'filled')
ylabel('mean gradient magnitude')

% -->> There is no clear linear dependence from the gradient to the number
%      of seeds. There is a trend that for more seeds, the gradients of
%      single seeds is smaller.


%% Where are the relevant seeds?
% Seed considered relevant are the seeds which have a gradient magnitude
% which is an outlier in the boxplots (figure 1).
% Points are outliers if they are greater than q3 + 1.5 × (q3 - q1) or less 
% than q1 - 1.5 × (q3 - q1), where q1 and q3 are the 25th and 75th 
% percentiles of the sample data, respectively.
% a) look at positioning in cartesian coordinate system
% b) look at position in spherical coordinate system
relXPositions = [];
relYPositions = [];
relZPositions = [];
dCoMOutliers  = [];
dCoMOther = [];
for i = 1:size(patientIds,2)
    % find outliers
    q = quantile(data{i}.V100_gradMagnitudes,[0.25 0.75]);
    upperWhiskerLimit = q(2) + 1.5 * (q(2) - q(1));
    data{i}.outlierIdxs = find(data{i}.V100_gradMagnitudes > upperWhiskerLimit);
    % b) reference
    for s = 1:size(data{i}.tplan,1)
        seedPosition = data{i}.tplan{s,2};
        if(ismember(s, data{i}.outlierIdxs))
            relXPositions = [relXPositions, (seedPosition(1) - data{i}.xExtendProstate(1))/(data{i}.xExtendProstate(2) - data{i}.xExtendProstate(1))];
            relYPositions = [relYPositions, (seedPosition(2) - data{i}.yExtendProstate(1))/(data{i}.yExtendProstate(2) - data{i}.yExtendProstate(1))];
            relZPositions = [relZPositions, (seedPosition(3) - data{i}.zExtendProstate(1))/(data{i}.zExtendProstate(2) - data{i}.zExtendProstate(1))];
            % b)
            dCoMOutliers = [dCoMOutliers, norm(seedPosition' - data{i}.centerOfMass)];
        else
            dCoMOther = [dCoMOther, norm(seedPosition' - data{i}.centerOfMass)];
        end
    end
end
ga = [repmat({'x'}, size(relXPositions,2),1); repmat({'y'}, size(relXPositions,2),1); repmat({'z'}, size(relXPositions,2),1)];
gb = [repmat({'other seeds'}, size(dCoMOther,2),1); repmat({'most relevant seeds'}, size(dCoMOutliers,2),1)];

figure
boxplot([relXPositions';relYPositions'; relZPositions'],ga);
ylabel('position in cartesian coordinates relative to prostate extend')
figure
boxplot([dCoMOther';dCoMOutliers'],gb);
ylabel('distance to center of mass')


% -->> a) The most relevant seeds are equally distributed in x, y, and z
%      relative to the prostate volume.
%      b) Similar in distance to center of mass

%% Multivariate linear regression
maxGradMag = [];
meanGradMag = [];
numSeeds = [];
prostateVol = [];
variableValue = [];
for i = 1:size(patientIds,2)
      meanGradMag = [meanGradMag; data{i}.mean_gradMagnitudes];
      maxGradMag  = [maxGradMag; data{i}.max_gradMagnitudes];
      numSeeds = [numSeeds; size(data{i}.tplan,1)];
      prostateVol = [prostateVol; data{i}.prostateVolume/1000];
      variableValue = [variableValue;data{i}.V100];
end
x = [numSeeds, prostateVol, variableValue];
model_mean = fitlm(x,meanGradMag)
figure
plot(model_mean)