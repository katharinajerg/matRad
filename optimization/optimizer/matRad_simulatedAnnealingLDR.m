function [x,fVal,exitflag,output] = matRad_simulatedAnnealingLDR(fHandle,x0,noSeeds,options)
%SIMULATEDANNEALINGLDR Summary of this function goes here
%   Detailed explanation goes here

startTime = tic;

switch nargin
    case 3
        options = [];
end

if isempty(options)
    options.tolFun          = 1e-3;
    options.maxIter         = 0;
    options.stallIterLimit  = 50;
    options.timeLimit       = 600;
    
    options.additionalSeeds = 0;
    options.additionalPositions = 0;
    options.maxTempIter = 500;
    options.linearParamTemperature = 10;
    options.speedParam = 0.6;

    options.importDefaultSeedPos = 1; %1 for seed position import. When importing seeds options.additialSeeds/Pos must be 0
end

if options.additionalSeeds<=1
    totSeeds = round(options.additionalSeeds*noSeeds)+noSeeds;
else
    totSeeds = options.additionalSeeds+noSeeds;
end

noPositions = length(x0);
if options.additionalPositions<=1
    totPositions = round(options.additionalPositions*noPositions)+noPositions;
else
    totPositions = noPositions+length(x0);
end

supp     = getSpacedSeedConfiguration(totPositions,totSeeds, options);
tempSupp = supp(1:noPositions);
fVal     = fHandle(tempSupp);

%% getting the initial temperature
initTemp = 0;
oldFVal  = fVal;
for i=1:options.maxTempIter
    supp = getNextSeedConfiguration(supp,totPositions);
    tempSupp = supp(1:noPositions);
    fVal = fHandle(tempSupp);
    initTemp = initTemp+abs(oldFVal-fVal);
end
initTemp = options.linearParamTemperature*initTemp/options.maxTempIter;

%% initialization of mandatory variables for SA
supp     = getSpacedSeedConfiguration(totPositions,totSeeds, options);
tempSupp = supp(1:noPositions);
fVal     = fHandle(tempSupp);

iter            = 0;
noStallIter     = 0;
funccount       = options.maxTempIter;
actFoundMinFVal = fVal;
actFoundMinSupp = supp;

while iter<options.maxIter && toc(startTime)<options.timeLimit && noStallIter<options.stallIterLimit
    iter=iter+1;
    actTemp  = initTemp/iter^options.speedParam;
    
    oldSupp  = supp;
    supp     = getNextSeedConfiguration(supp,totPositions);
    tempSupp = supp(1:noPositions);
    tempFVal = fHandle(tempSupp);
    funccount = funccount+1;
    dFVal    = tempFVal-fVal;
    
    if dFVal<=0
        fVal = tempFVal;
        
        if fVal<actFoundMinFVal
            actFoundMinFVal = fVal;
            actFoundMinSupp = supp;
        end
    else
        valueToCheck      = rand();
        valueToBeAccepted = exp(-dFVal/actTemp);
        
        if valueToCheck<valueToBeAccepted
            fVal = tempFVal;
        else
            supp = oldSupp;
        end
    end

    
    if (dFVal<options.tolFun)
        noStallIter=noStallIter+1;
    else
        noStallIter=0;
    end
     
end
    
lastToc = toc(startTime);
supp = actFoundMinSupp;
fVal = actFoundMinFVal;


if iter==options.maxIter
    exitflag=-2;
end
if noStallIter==options.stallIterLimit
    exitflag=0;
end
if lastToc>options.timeLimit
    exitflag=-1;
end


fprintf('* Q(%i) = %4.2f   (elapsed time t = %4.2f) \n',iter,fVal,lastToc);
fprintf('* ')
% for i=1:length(supp)
%     fprintf('%i\t',supp(i))
%     if mod(i,10)==0
%         fprintf('\n* ')
%     end
% end
fprintf('\n')

x = supp(1:noPositions);

output.iterations = iter;
output.funccount  = funccount;
output.totaltime  = lastToc;
output.algorithm  = 'Simulated Annealing';
output.exitflag    = exitflag;


end


function supp = getSpacedSeedConfiguration(totPositions,totSeeds, options)

    if options.importDefaultSeedPos
        load('supp.mat', 'supp');
        supp = double(supp);
    else
        distSeed = floor(totPositions/(totSeeds+2));
        ind = distSeed:distSeed:totSeeds*distSeed;
        supp = zeros(totPositions,1);
        supp(ind) = 1;
    end
end


function supp = getNextSeedConfiguration(supp,maxPos)
    ind = find(supp);
    indSupp = randi(numel(ind));
    dir     = 2*(randi(2)-1)-1;
    
    indexToAdd = ind(indSupp);
    validConfig = 0;

    while validConfig==0
        indexToAdd = indexToAdd+dir;
        if indexToAdd<1
            indexToAdd = maxPos;
        end
        if indexToAdd>maxPos
            indexToAdd = 1;
        end
    
        if ~any(ind==indexToAdd)
            validConfig=1;
        end
    end

    supp(ind(indSupp)) = 0;
    supp(indexToAdd) = 1;

end



