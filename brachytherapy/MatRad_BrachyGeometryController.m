classdef MatRad_BrachyGeometryController
    %MATRAD_BRACHYGEOMETRYCONTROLLER Handles the calculation from  a Set of
    %support point sets to its dwell points.
    %   For each support point set it gets dwell points by
    %   using the MatRad_BrachyNeedle class.
    
    properties
        supportPointSet     % Support point set for all needles
        numberOfNeedles     % Total number of Needles
        needleSet           % Set of all needles
    end
    
    methods
        function obj = MatRad_BrachyGeometryController(mySuppPointSet, myNumberOfNeedles)
            %MATRAD_BRACHYGEOMETRYCONTROLLER Construct an instance of this class
            %    Setting base parameters
            obj.supportPointSet = mySuppPointSet;
            obj.numberOfNeedles = myNumberOfNeedles;
            obj.needleSet = {};
          
        end
        
        
        function [obj] = calcNeedles(obj)
            %Each support point set is handled by an instance  of
            %MatRad_BrachyNeedle.
            myNeedles  = {obj.numberOfNeedles};
            parfor n = 1:obj.numberOfNeedles
                n;
                needle = MatRad_BrachyNeedle(obj.supportPointSet(n), 15, 3, 5);
                needle = needle.calc();
                

                myNeedles{n} = needle;
                
            end
            obj.needleSet = myNeedles;
            
        end
        
    end
end

