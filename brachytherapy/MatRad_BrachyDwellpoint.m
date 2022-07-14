classdef MatRad_BrachyDwellpoint
    %MATRAD_BRACHYDWELLPOINTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        positionX
        positionY
        positionZ
        orientationX
        orientationY
        orientationZ
        length
    end
    
    methods
        function obj = MatRad_BrachyDwellpoint(myPosition, myOrientation)
            %MATRAD_BRACHYDWELLPOINTS
            obj.positionX = myPosition(1);
            obj.positionY = myPosition(2);
            obj.positionZ = myPosition(3);

            myOrientation = myOrientation/norm(myOrientation);
            obj.orientationX = myOrientation(1);
            obj.orientationY = myOrientation(2);
            obj.orientationZ = myOrientation(3);
            obj.length = 2;
        end
        
        function obj = setLength(obj, newLength)
            %Set Length
            %   Input: length
            obj.length = newLength;
        end

        function length = getLength(obj)
            %Get Length
            %   output: length
            length = obj.length;
        end

        function obj = setPosition(obj, newPosition)
            %Set position
            %   input: [x y z]
            obj.positionX = newPosition(1);
            obj.positionY = newPosition(2);
            obj.positionZ = newPosition(3);
        end

        function position = getPosition(obj)
            %Get position
            %   output: [x y z]
            position = [obj.positionX obj.positionY obj.positionZ];
        end

        function obj = setOrientation(obj, newOrientation)
            %Set nomalitzed orientation
            %   input: [x y z]
            newOrientation = newOrientation/norm(newOrientation);
            obj.positionX = newOrientation(1);
            obj.positionY = newOrientation(2);
            obj.positionZ = newOrientation(3);
        end

        function orientation = getOrientation(obj)
            %%Get orientation
            %   output: [x y z]
            orientation = [obj.orientationX obj.orientationY obj.orientationZ];
        end
    end
end

