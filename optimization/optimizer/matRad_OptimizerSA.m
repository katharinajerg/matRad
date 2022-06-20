classdef matRad_OptimizerSA < matRad_Optimizer
% matRad_OptimizerSA implements the interface for simulated annealing in
% LDR
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
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        options
        wResult
        resultInfo
        env
    end
    
    properties (Access = private)
        allObjectiveFunctionValues
        axesHandle
        plotHandle
        abortRequested
        plotFailed
    end
    
    methods
        function obj = matRad_OptimizerSA
            %matRad_OptimizerSA
            %   Construct an instance of the simulated annealing optimizer
            
            matRad_cfg = MatRad_Config.instance();
            
            obj.wResult = [];
            obj.resultInfo = [];
            obj.axesHandle = [];
            obj.allObjectiveFunctionValues = [];
            obj.abortRequested = false;
                     

        end
        
        function obj = optimize(obj,w0,optiProb,dij,cst)
            matRad_cfg = MatRad_Config.instance();
            
            % set optimization options            
            saStruct = struct();           

            % constraint bounds;
            [saStruct.cl,saStruct.cu] = optiProb.matRad_getConstraintBounds(cst);
            
            % set callback functions
            fHandle = @(x) optiProb.matRad_objectiveFunction(x,dij,cst);
%             fHandle       = @(x) optiProb.matRad_constraintFunctions(x,dij,cst);
                       

            % define number of seeds
            noSeeds = 30;%numel(w0);

            % set Callback
            qCallbackSet = false;
            if ~isdeployed % only if _not_ running as standalone                              
                try
                    % get handle to Matlab command window
                    mde         = com.mathworks.mde.desk.MLDesktop.getInstance;
                    cw          = mde.getClient('Command Window');
                    xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);
                    h_cw        = handle(xCmdWndView,'CallbackProperties');
                    
                    % set Key Pressed Callback of Matlab command window
                    set(h_cw, 'KeyPressedCallback', @(h,event) obj.abortCallbackKey(h,event));
                    fprintf('Press q to terminate the optimization...\n');
                    qCallbackSet = true;
                catch
                    matRad_cfg.dispInfo('Manual termination with q not possible due to failing callback setup.\n');
                end                
            end

            
            % run Simulated Annealing
            try
                [obj.wResult, ~ , exitflag, obj.resultInfo] = matRad_simulatedAnnealingLDR(fHandle,w0, noSeeds);
            catch ME
                errorString = [ME.message '\nThis error was thrown by simulated Annealing LDR.'];
                matRad_cfg.dispError(errorString);
            end
            
            % unset Key Pressed Callback of Matlab command window
            if qCallbackSet
                set(h_cw, 'KeyPressedCallback',' ');
            end
            
            obj.abortRequested = false;
            % Empty the array of stored function values
            obj.allObjectiveFunctionValues = [];
        end

               
        function [statusmsg,statusflag] = GetStatus(exitflag)
            try
                switch exitflag
                    case -2
                        statusmsg = 'maxIter';
                    case -1
                        statusmsg = 'time limit';
                    case 0
                        statusmsg = 'solved';
                    otherwise
                        statusmsg = 'SA returned unknown status';
                end
                
                if exitflag == 0
                    statusflag = 0;
                else 
                    statusflag = 1;
                end
                              
            catch
                statusmsg = 'No Last SA Status Available!';
                statusflag = -1;
            end
        end
        
        function plotFunction(obj)
            % plot objective function output
            y = obj.allObjectiveFunctionValues;
            x = 1:numel(y);
            
            if isempty(obj.axesHandle) || ~isgraphics(obj.axesHandle,'axes')
                %Create new Fiure and store axes handle
                hFig = figure('Name','Progress of SA Optimization','NumberTitle','off','Color',[.5 .5 .5]);
                hAx = axes(hFig);
                hold(hAx,'on');
                grid(hAx,'on');
                grid(hAx,'minor');
                set(hAx,'YScale','log');
                 
                %Add a Stop button with callback to change abort flag
                c = uicontrol;
                cPos = get(c,'Position');
                cPos(1) = 5;
                cPos(2) = 5;
                set(c,  'String','Stop',...
                        'Position',cPos,...
                        'Callback',@(~,~) abortCallbackButton(obj));                
                
                %Set up the axes scaling & labels
                defaultFontSize = 14;
                set(hAx,'YScale','log');
                title(hAx,'Progress of Optimization','LineWidth',defaultFontSize);
                xlabel(hAx,'# iterations','Fontsize',defaultFontSize),ylabel(hAx,'objective function value','Fontsize',defaultFontSize);
                
                %Create plot handle and link to data for faster update
                hPlot = plot(hAx,x,y,'xb','LineWidth',1.5,'XDataSource','x','YDataSource','y');
                obj.plotHandle = hPlot;
                obj.axesHandle = hAx;
                                
            else %Figure already exists, retreive from axes handle
                hFig = get(obj.axesHandle,'Parent');
                hAx = obj.axesHandle;
                hPlot = obj.plotHandle;
            end
          
            drawnow;
            
            % ensure to bring optimization window to front
            figure(hFig);
        end
        
        function abortCallbackKey(obj,~,KeyEvent)
            % check if user pressed q
            if  get(KeyEvent,'keyCode') == 81
                obj.abortRequested = true;
            end
        end
        
        function abortCallbackButton(obj,~,~,~)
            obj.abortRequested = true;
        end        
    end
    
    methods (Static)
        function available = IsAvailable()
            available = matRad_checkMexFileExists('ipopt');                   
        end
    end
end
