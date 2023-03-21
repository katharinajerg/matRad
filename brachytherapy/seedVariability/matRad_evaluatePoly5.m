function funValues = matRad_evaluatePoly5(coeff, xValues)
%MATRAD_EVALUATEPOLY5 evaluates a polynomial of degree 5:
% y = a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f 
%
% input
%   coeff:  coefficient vector [a,b,c,d,e,f] 
%   x:      vector of x values to be evaluated
%	
% output
%   y:      vector of function values    

    funValues = coeff(1).* xValues.^5 + coeff(2).*xValues.^4 + coeff(3).*xValues.^3 ...
        + coeff(4).*xValues.^2 + coeff(5).*xValues + coeff(6);
end


