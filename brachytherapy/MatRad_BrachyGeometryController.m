classdef MatRad_BrachyGeometryController
    %MATRAD_BRACHYGEOMETRYCONTROLLER Handles the calculation from  a Set of
    %support point sets to its dwell points.
    %   For each support point set it gets dwell points by
    %   using the MatRad_BrachyNeedle class.
    
    properties
        supportPointSet     % Support point set for all needles
        numberOfNeedles     % Total number of Needles
        needleSet           % Set of all needles
        seedDistance        % Distance in mm between seeds
        numSeedsPerNeedle   % Number of seeds per needle
    end
    
    methods
        function obj = MatRad_BrachyGeometryController(mySuppPointSet, myNumberOfNeedles, mySeedDistance, myNumSeeds)
            %MATRAD_BRACHYGEOMETRYCONTROLLER Construct an instance of this class
            %    Setting base parameters
            obj.supportPointSet = mySuppPointSet;
            obj.numberOfNeedles = myNumberOfNeedles;
            obj.needleSet = {};
            obj.seedDistance = mySeedDistance;
            obj.numSeedsPerNeedle = myNumSeeds;

          
        end
        
        
        function [obj] = calcNeedles(obj)
            %Each support point set is handled by an instance  of
            %MatRad_BrachyNeedle.
            myNeedles  = {obj.numberOfNeedles};
            for n = 1:obj.numberOfNeedles
                needle = MatRad_BrachyNeedle(obj.supportPointSet(n), 200, obj.seedDistance, obj.numSeedsPerNeedle);
                needle = needle.calc();
                
                myNeedles{n} = needle;
                
            end
            obj.needleSet = myNeedles;
            
        end
        
    end
end

