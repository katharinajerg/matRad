function [OrientationMatrix,OrientationVector] = matRad_getOrientationMatrix(seedPoints,dosePoints)
%MATRAD_GETORIENTATIONMATRIX gets (seedpoint x dosepoint) matrix of relative
% distances

% input
%   seedPoints:     struct with fields x,y,z
%   dosePoints:     struct with fields x,y,z
%
% output
%   orientation matrix:    rows: index of dosepoint 
%                       columns: index of seedpoint
%                       entry: orientation of seedpoint
%                       
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
OrientationMatrix.x_orientation = ones(1,length(dosePoints.x))'*seedPoints.x_orientation;
OrientationMatrix.y_orientation = ones(1,length(dosePoints.y))'*seedPoints.y_orientation;
OrientationMatrix.z_orientation = ones(1,length(dosePoints.z))'*seedPoints.z_orientation;

if nargout == 2
OrientationVector = reshape(OrientationMatrix.dist,[],1);
end

end

