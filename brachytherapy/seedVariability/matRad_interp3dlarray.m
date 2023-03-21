function Y = matRad_interp3dlarray(xi,yi,zi,I,xq,yq,zq)
%MATRAD_INTERP3DLARRAY Trilinear interpolation for data in dlarray
% Attention: the function is a function I(y,x)!!
%
% call
%   y = matRad_interp3dlarray(xi,yi,zi,I,xq,yq,zq)
%
% input
%   xi:         x coordinates of the sample points 
%   yi:         y coordinates of the sample points 
%   zi:         z coordinates of the sample points 
%   I:          corresponding function values in a dlarray  
%   xq:         x coordinates of the query points
%   yq:         y coordinates of the query points
%   zq:         z coordinates of the query points
%
% output
%   Y:          function values at the query points in dlarray
%
% references: 
% interpolation weights from:
% https://en.wikipedia.org/wiki/Trilinear_interpolation
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
    Y = dlarray(zeros(numel(yq),numel(xq),numel(zq)));

    for k = 1:numel(zq)
        z = zq(k);
        for j=1:numel(yq)
            y = yq(j);
            for i=1:numel(xq)
                x = xq(i);
          % Any values out of acceptable range
              x(x < xi(1)) = xi(1);
              x(x > xi(end) - 0.001) = xi(end) - 0.001;
              y(y < yi(1)) = yi(1);
              y(y > yi(end) - 0.001) = yi(end) - 0.001;
              z(z < zi(1)) = zi(1);
              z(z > zi(end) - 0.001) = zi(end) - 0.001;
          % get indicies for x
              x1 = 1;
              for n = 2:numel(xi)
                  if xi(n)<x
                      x1 = n;
                  else 
                      break;           
                  end
              end 
              x2 = x1 + 1;

              y1 = 1;
              for n = 2:numel(yi)
                  if yi(n)<y
                      y1 = n;
                  else 
                      break;           
                  end
              end    
              y2 = y1 + 1;

              z1 = 1;
              for n = 2:numel(zi)
                  if zi(n)<z
                      z1 = n;
                  else 
                      break;           
                  end
              end
              z2 = z1 + 1;

          % 8 Neighboring Pixels
              NP1 = I(y1,x1,z1);
              NP2 = I(y1,x2,z1);
              NP3 = I(y2,x1,z1); 
              NP4 = I(y2,x2,z1);
              NP5 = I(y1,x1,z2);
              NP6 = I(y1,x2,z2);
              NP7 = I(y2,x1,z2);
              NP8 = I(y2,x2,z2);

          % Interpolation weights 
              xd = (x-xi(x1))/(xi(x2)-xi(x1));
              yd = (y-yi(y1))/(yi(y2)-yi(y1));
              zd = (z-zi(z1))/(zi(z2)-zi(z1));

              c11 = NP1*(1-xd) + NP2*xd;
              c12 = NP5*(1-xd) + NP6*xd;
              c21 = NP3*(1-xd) + NP4*xd;
              c22 = NP7*(1-xd) + NP8*xd;

              c1 = c11*(1-yd) + c21*yd;
              c2 = c12*(1-yd) + c22*yd;

          % Function value
              Y(j,i,k) = c1*(1-zd) + c2*zd;
            end
        end
    end
end
