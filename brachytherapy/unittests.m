classdef unittests < matlab.unittest.TestCase
    
    methods(Test)
        % Test methods

        function testDwellpoint(testCase)
            myDwellpoint = MatRad_BrachyDwellpoint([1 2 -1], [4 5 6]);
            testCase.verifyEqual(myDwellpoint.positionX, 1);
            testCase.verifyEqual(myDwellpoint.positionY, 2);
            testCase.verifyEqual(myDwellpoint.positionZ, -1)
            testCase.verifyEqual(myDwellpoint.getPosition, [1 2 -1]);
            testCase.verifyEqual(size(myDwellpoint.getPosition), [1 3]);
            n = norm([4 5 6]);

            testCase.verifyEqual(myDwellpoint.orientationX, 4/n);
            testCase.verifyEqual(myDwellpoint.orientationY, 5/n);
            testCase.verifyEqual(myDwellpoint.orientationZ, 6/n);
            testCase.verifyEqual(myDwellpoint.getOrientation, [4 5 6]/n);
            testCase.verifyEqual(size(myDwellpoint.getOrientation), [1 3]);

            testCase.verifyEqual(myDwellpoint.length, 2);
        end


        function testNeedleCalc(testCase)
            supportPoints ={[
                                 0 2 0; 1 4 0; 2.5 6 0; 3.6 7 0; 5 8 0; 7 8.75 0; 8.1 8.89 0; 10 9 0; 12 9.9 0;
                            ]};
                           
            myNeedle =  MatRad_BrachyNeedle(supportPoints, 15, 3, 5);
            testCase.verifyEqual(myNeedle.length, 15);
            testCase.verifyEqual(myNeedle.distance, 3);
            testCase.verifyEqual(myNeedle.numberOfSeeds, 5);
            testCase.verifyEqual(size(myNeedle.supportPoints{1}), [9 3]);

            
            myNeedle = myNeedle.checkEnvironment;
            
            testCase.verifyEqual(myNeedle.environment, 'Plane');
            testCase.verifyEqual(size(myNeedle.dwellPoints), [0 0]);
            testCase.verifyEqual(size(myNeedle.dwellPointList), [0 0]);

            myNeedle = myNeedle.calcEnvironment;
            testCase.verifyEqual(size(myNeedle.dwellPoints), [5 3]);
            testCase.verifyEqual(size(myNeedle.dwellPointList), [1 5]);
            
        end
        
        function testController(testCase)

            load data/suppPoints1.mat
            myController = MatRad_BrachyGeometryController(supportPoints, 4);

            testCase.verifyEqual(myController.numberOfNeedles, 4);
            testCase.verifyEqual(size(myController.supportPointSet), [4 1]);
            testCase.verifyEqual(size(myController.needleSet), [0 0]);

            myController = myController.calcNeedles;
            testCase.verifyEqual(size(myController.needleSet), [1 4]);
            testCase.assertClass(myController.needleSet{1}, 'MatRad_BrachyNeedle');
            testCase.assertClass(myController.needleSet{2}, 'MatRad_BrachyNeedle');
            testCase.assertClass(myController.needleSet{3}, 'MatRad_BrachyNeedle');
            testCase.assertClass(myController.needleSet{4}, 'MatRad_BrachyNeedle');

            

        end
    end
    
end