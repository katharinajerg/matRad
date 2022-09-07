function [ThetaMatrix,ThetaVector] = matRad_getThetaMatrix(OrientationMatrix,DistanceMatrix)
% getThetaMatrix gets (seed x dosepoint) matrix of relative polar angles
%
% call
%   [ThetaMatrix,ThetaVector] = matRad_getThetaMatrix(OrientationMatrix,...
%       DistanceMatrix)
%   normally called within matRad_calcBrachyDose
%   !!getDistanceMatrix needs to be called first!!
%   !!getOrientationMatrix needs to be called first!!
%
% input
%   DistanceMatrix:     [dosePoint x seedPoint] struct with fields 'x','y',
%                       'z' and total distance 'dist'
%   OrientationMatrix:  [dosePoint x seedPoint] struct with fields 'x','y',
%                       'z', orientation for each seed point and its dose
%                       points are equaly
%
% output
%   angle matrix:       rows: index of dosepoint 
%                       columns: index of deedpoint
%                       entry: polar angles between seedpoints and  
%                       dosepoint in degrees
%   angle vector:       column vector of angle matrix entries
%
% comment:
%   The shape of the Theta matrix will be consistent with the shape of 
%   input fields. 
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


DistanceMatrix.dist(DistanceMatrix.dist == 0) = 1; %Avoid deviding by zero
ThetaMatrix = acosd((OrientationMatrix.x .*DistanceMatrix.x ...
                        + OrientationMatrix.y .*DistanceMatrix.y ...
                        + OrientationMatrix.z .*DistanceMatrix.z)./DistanceMatrix.dist);  
if nargout == 2
    ThetaVector = reshape(ThetaMatrix,[],1);
end

end

