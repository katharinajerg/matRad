function [xSorted] = matRad_sortdlarray(x)
%MATRAD_SORTDLARRAY sort values inside a dlarray in ascending order
% 
% input
%   x:      vector of numbers as a dlarray
%	
% output
%   xSorted: dlarray with sorted values of input vector x
%
% references
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

    xSorted = x;
    for i = 1:numel(x)
        [minValue, ind] = min(x);
        xSorted(i) = minValue;
        x(ind) = [];   
    end
end

