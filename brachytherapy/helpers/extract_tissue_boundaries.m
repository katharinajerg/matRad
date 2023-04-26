
function [boundaries] = extract_tissue_boundaries(ct,cst)
%EXTRACT_TISSUE_BOUNDARIES Summary of this function goes here
%   Detailed explanation goes here
% find all target voxels from cst cell array
boundaries = cell(size(cst,1),2);
for t = 1:size(cst,1)
  
    V = vertcat(cst{t,4}{:});
    
    % Remove double voxels
    V = unique(V);
    % generate voi cube for targets
    voiTarget    = zeros(ct.cubeDim);
    voiTarget(V) = 1;
        
  
    % Convert linear indices to 3D voxel coordinates
    [coordsY_vox, coordsX_vox, coordsZ_vox] = ind2sub(ct.cubeDim,V);
    
    % take only voxels inside patient
    V = [cst{:,4}];
    V = unique(vertcat(V{:}));
    
    % ignore densities outside of contours
    eraseCtDensMask = ones(prod(ct.cubeDim),1);
    eraseCtDensMask(V) = 0;
    for i = 1:ct.numOfCtScen
        ct.cube{i}(eraseCtDensMask == 1) = 0;
    end
    
    % save VOI coordinates
    %translate to geometric coordinates and save in stf  
    TargX = ct.x(coordsX_vox);
    TargY = ct.y(coordsY_vox);
    TargZ = ct.z(coordsZ_vox);
    
    P = [TargX',TargY',TargZ'];
    k = boundary(P,1);

    boundaries{t,1} = cst{t,2};
    boundaries{t,2} = k;
    boundaries{t,3} = P;
    
end

end

