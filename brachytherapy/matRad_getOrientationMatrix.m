function [OrientationMatrix,OrientationVector] = matRad_getOrientationMatrix(seedPoints,dosePoints)
%MATRAD_GETORIENTATIONMATRIX Summary of this function goes here
%   Detailed explanation goes here

OrientationMatrix.x_orientation = ones(1,length(dosePoints.x))'*seedPoints.x_orientation;
OrientationMatrix.y_orientation = ones(1,length(dosePoints.y))'*seedPoints.y_orientation;
OrientationMatrix.z_orientation = ones(1,length(dosePoints.z))'*seedPoints.z_orientation;

if nargout == 2
OrientationVector = reshape(OrientationMatrix.dist,[],1);
end

end

