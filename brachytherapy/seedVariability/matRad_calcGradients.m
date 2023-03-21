function [fval,gradval] = matRad_calcGradients(x)
    %MATRAD_calcGradient is a matRad function which calculates the
    % gradients of the quality measures for each seed position. 
    %
    % In order to automatically differentiate a function the input
    % needs to be a dlarray. The only purpose of this function is that, dlfeval 
    % of the function calculating the gradients (matRad_automaticDifferentiation)
    % must be called within a function. (See:
    % https://de.mathworks.com/help/deeplearning/ug/include-automatic-differentiation.html)
    %
    % input
    % x:        Vector of seed coordinates [x1,y1,z1, ..., xn, yn, zn] for 
    %           n seeds
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

    x0 = dlarray(x);
    [fval,gradval] = dlfeval(@matRad_automaticDifferentiation,x0);
end



