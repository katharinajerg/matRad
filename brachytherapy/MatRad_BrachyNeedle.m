classdef MatRad_BrachyNeedle
    %
    %   This class is used to calulate the geometrie 
    %   of needles needed for seed geometrie in Brachy.
    %
    %   In this approach the support points are expected
    %   to be in line or fitting in one plane. 
    %   So in the first step it fits a plane that all points are
    %   in the plane. With the plane and its norm vector there is a now 
    %   new coordinate system. the origin is set to the first support
    %   point. All support points are in the x1,x2 plane so we can solve 
    %   curve fitting and arc length problems in 2D.
    %   Resulting dwell points are stored in the 'dwellPointList' variable.
    %   It contains instances of the MatRad_BrachyDwellpoint class 
    
    properties
        supportPoints           % Points where the needle should go through
        environment             % 'Line' or 'Plane', if support points are in a Line or not.

        dwellPoints             % resulting dwell points coordinates
        dwellPointList          % resulting dwell point List, containing position and orientation

        planePoints             % Points in plane coordinats
        dwellPointsPlane        % resulting dwell point in plane coordinates 
        dwellPointsPlaneList    % resulting dwell point List, containing position and orientation in plane coordinates 

        length                  % length of the needle
        distance                % distance between the Seeds  
        numberOfSeeds           % amount of Seed pro Needle

        plane                   % 3 vectors to describe teh plane [ Stützvektor, Spannvektor1, Spannvktor2] 
        T                       % change-of-basis Matrix plane to standard basis
        T1                      % change-of-basis Matrix plane to standard basis
        spline                  % piecewise polynomial fitting planepoints
    end
    
    methods
        function obj = MatRad_BrachyNeedle(suppSet,myLength, myDistance, myNumOfSeeds)
            %MatRad_Needle Construct an instance of this MatRad_BrachyNeedle class
            %   Setting base parameters

            obj.supportPoints = suppSet;
            obj.length = myLength;
            obj.distance = myDistance;
            obj.numberOfSeeds = myNumOfSeeds;
            obj.T = ones(3,3);
        end


        function obj = calc(obj)
            %Calculates dwell points for the needle with its given support
            %points. By calling member functions.

            obj = obj.checkEnvironment();
            obj = obj.calcEnvironment();
            
        end

        function obj = checkEnvironment(obj)
            %Verify if the support points are in a plane or line
            %and prepares necessary variables.
            %The result is stored in the class variable 'environment'
            
            thirdPointIndex = 3;
            SuppPointsDim = (size(obj.supportPoints{1}));
            anzSuppPoints = SuppPointsDim(1);
            e1 = obj.supportPoints{1}(1, :);
            e2 = obj.supportPoints{1}(2, :);
            e3 = obj.supportPoints{1}(thirdPointIndex, :);
            v1 = e2 -e1;
            v2 = e3 -e1;
            [v1, v2];

            while matRad_linearDependent(v1,v2) && (thirdPointIndex < anzSuppPoints)
                %'Check next Point';
                thirdPointIndex  = thirdPointIndex+1;
                e3 = obj.supportPoints{1}(thirdPointIndex, :);
                v2 = e3 -e1;
            end

            if matRad_linearDependent(v1,v2)
                %'line';
                obj.environment = 'Line';
            else
                %'plane';
                v1 = v1/norm(v1);
                obj.environment = 'Plane';
                n = cross(v1, v2);
                n = n/norm(n);
                v2 = cross(n, v1);
                v2 = v2/norm(v2);
                myPlane = {e1 v1 v2};
                obj.T(:,1) = v1;
                obj.T(:,2) = v2;
                obj.T(:,3) = n;
                obj.T;
                obj.plane = myPlane;
                obj.T1 = inv(obj.T);
                obj.T1;
            end
        end

        function obj = calcEnvironment(obj)
            % for its enviroment it calls the corresponding member function
            % to calculate the dwellpoints and its orientations.
            
            if strcmp(obj.environment , 'Line')
                %gerade
                'Line';
                 obj = obj.calcLine();
            elseif strcmp(obj.environment,'Plane')
                %ebene
                obj = obj.calcPlanePoints();
                obj = obj.calcDistances();
                obj = obj.calcDwellPoints();
            end
        end

        function obj = calcPlanePoints(obj)
            %specifies plane coordinates for all support points.
            % resulting points should be all in a x1-x2 plane.
            anzPoints = size(obj.supportPoints{1});
            obj.planePoints = zeros([anzPoints(1) 2]);

            for i= 1:anzPoints(1)
                planePoint = obj.T1 *transpose(obj.supportPoints{1}(i, :)-obj.plane{1});
                obj.planePoints(i, :) = planePoint(1:2);
                if planePoint(3)~=0
                   %'not plane'
                    assert(planePoint(3)==0, "Supportpunkte liegen nicht auf einer Ebene! Geben sie andere Punkte an.");
                end
            end
        end

        function obj = calcDistances(obj)
            %determines dwellpoint with the right distance by considering
            %the arc length

            obj.spline = csapi(obj.planePoints(:, 1), obj.planePoints(:, 2));
            
            [breaks,coefs,l,k,d] = unmkpp(obj.spline);
            f = @(x)ppval(obj.spline, x);
            % make the polynomial that describes the derivative
            diff = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
            diffn = @(x)ppval(diff, x);
            diffn2 = @(x)sqrt(1+diffn(x)^2);
            blfn = @(b)integral(diffn2,0,b, "ArrayValued",true);
          
            anzPoints = size(obj.planePoints);
            maxSeets = obj.numberOfSeeds;

            mydwellPointsdistancesPlane = zeros([anzPoints(1) 2]);
            seedPointsPlane = zeros([maxSeets 2]);
            obj.dwellPointsPlane = zeros([maxSeets 2]);

            i = 1;
            v = 1;
            obj.dwellPointsPlaneList = struct([]);

            %Ausgangspunkt
            fdistance = blfn(0);
            f0 = f(0);
            seedPointsPlane(v, : ) =  [0, f0];
            obj.dwellPointsPlane(v, : ) =  [0,f0];
            obj.dwellPointsPlaneList{v} = {[0,f0], diffn(0)};
            for x= 0.2:0.2:obj.planePoints(end, 1)
                if(v >= maxSeets)
                    break
                end
                fdistance = blfn(x);
                fx=  f(x);
                mydwellPointsdistancesPlane(i, :) = [x, fdistance(1)];
                if(fdistance(1)-seedPointsPlane(v,2) > obj.distance)
                    v = v+1;
                    seedPointsPlane(v, : ) =  [x, fdistance(1)];
                    obj.dwellPointsPlane(v, : ) = [x, fx(1)];
                    obj.dwellPointsPlaneList{v} = {[x, fx(1)], diffn(x)};
                end
                i = i+1;
            end
            if(v < maxSeets)
                warning('Needle is too short to fit all wanted dwell points. Check input data.')
            end
        end

        function obj = calcDwellPoints(obj)
            %back transformation to standard basis
            
            obj.planePoints;
            
            anzPoints = size(obj.dwellPointsPlane);
            obj.dwellPoints = zeros([anzPoints(1) 3]);
            for i= 1:anzPoints(1)
                p = [obj.dwellPointsPlane(i, :) 0];
                dwellPoint = obj.T * transpose(p) + transpose(obj.plane{1});
                m = (obj.dwellPointsPlaneList{i}(1,2));
                v1 = obj.T * transpose([1 m{1} 0]);
                v1 = v1/norm(v1);
                obj.dwellPoints(i, :) = dwellPoint;
                obj.dwellPointList{i} = MatRad_BrachyDwellpoint(dwellPoint, v1);
            end
        end

        function obj = checkT(obj)           
            anzPoints = size(obj.planePoints);
            obj.dwellPoints = zeros([anzPoints(1) 3]);
            for i= 1:anzPoints(1)
                i
                p = [obj.planePoints(i, :) 0];
                dwellPoint = obj.T * transpose(p) + transpose(obj.plane{1});
                sp = obj.supportPoints{1}(i, :);
                assert(sp(1) == dwellPoint(1), strcat('T nicht korrekt! ', int2str(sp(1)),'~= ' , int2str(dwellPoint(1))))
                assert(sp(2) == dwellPoint(2),  strcat('T nicht korrekt! ', int2str(sp(2)),'~= ' , int2str(dwellPoint(2))))
                assert(sp(3) == dwellPoint(3),  strcat('T nicht korrekt! ', int2str(sp(3)),'~= ' , int2str(dwellPoint(3))))
            end
        end

        function obj = calcLine(obj)
            %determines dwellpoint if support points are in line
            obj.dwellPointList = struct([]);
            e1 = obj.supportPoints{1}(1, :);
            e2 = obj.supportPoints{1}(2, :);
            v1 = e2 -e1; 
            v1 = v1/norm(v1);
            for seed = 0:obj.numberOfSeeds
                dwellPoint = obj.supportPoints{1}(1, :) +  seed*obj.distance*v1;
                obj.dwellPoints(seed+1, :) = transpose(dwellPoint);            
                obj.dwellPointList{seed+1} = MatRad_BrachyDwellpoint(dwellPoint, v1);
            end
        end

    end
end

