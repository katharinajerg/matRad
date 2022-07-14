function bool = matRad_linearDependent(v1,v2)
%LINEARDEPENDENT checking two 2-d voctor for linear dependence
%   returns True if they are linear dependence
%
%   called by MatRad_BrachyNeedle.checkEnvironment()

    bool = v1(1)/v2(1) == v1(2)/v2(2) &&  v1(2)/v2(2) == v1(3)/v2(3);

end

