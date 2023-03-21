function [Y] = matRad_interp2dlarray(xi,yi,I,xq,yq)
%matRad_interp2dlarray returns interpolated values of a function of two 
% variables at specific query points using linear interpolation.  
% This function takes dlarrays as input, which the default MATLAB
% implementation (interp2) does not. 
% Attention: the function is a function I(y,x)!!
% 
% input
%   xi:     x coordinates of the sample points
%   yi:     y coordinates of the sample points
%   I:      function values at each sample point
%   xq:     x coordinates of the query points
%   yq:     y coordinates of the query points
%	
% output
%   Y: interpolated data   
%
%   References
%     -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    h = numel(yi);
    w = numel(xi);
    Y = zeros(numel(yq),numel(xq));
    hs = ((yq(2)-yq(1))/(yi(2)-yi(1)));
    ws = ((xq(2)-xq(1))/(xi(2)-xi(1)));
    ratio_x = 1/hs;
    ratio_y = 1/ws;
    for i=1:numel(yq)
      y = yq(i);
        for j=1:numel(xq)
           x = xq(j);
      % Any values out of acceptable range
          x(x < xi(1)) = xi(1);
          x(x > xi(end) - 0.001) = xi(end) - 0.001;
          x1 = floor(x);
          x2 = x1 + 1;
          y(y < yi(1)) = yi(1);
          y(y > yi(end) - 0.001) = yi(end) - 0.001;
          y1 = floor(y);
          y2 = y1 + 1;

      % 4 Neighboring Pixels
          NP1 = I(y1,x1);
          NP2 = I(y1,x2);
          NP3 = I(y2,x1); 
          NP4 = I(y2,x2);

      % 4 Pixels Weights
          PW1 = (y2-y)*(x2-x);
          PW2 = (y2-y)*(x-x1);
          PW3 = (x2-x)*(y-y1);
          PW4 = (y-y1)*(x-x1);
          
      % Function value    
          Y(i,j) = PW1 * NP1 + PW2 * NP2 + PW3 * NP3 + PW4 * NP4;
        end
    end
end

