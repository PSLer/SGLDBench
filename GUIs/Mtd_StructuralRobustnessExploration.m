classdef Mtd_StructuralRobustnessExploration < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        EvaluateDesignsStiffnesswrtChangedLoadingDirectionsPanel  matlab.ui.container.Panel
        ResultDisplayPanel              matlab.ui.container.Panel
        DesignCompliancewrtChangedLoadingDirectionsEditField  matlab.ui.control.NumericEditField
        DesignCompliancewrtChangedLoadingDirectionsEditFieldLabel  matlab.ui.control.Label
        DesignVolumeFractionEditField   matlab.ui.control.NumericEditField
        DesignVolumeFractionEditFieldLabel  matlab.ui.control.Label
        DesignComplianceEditField       matlab.ui.control.NumericEditField
        DesignComplianceEditFieldLabel  matlab.ui.control.Label
        RotateOriginalLoadingDirectionsviaEulersAnglesPanel  matlab.ui.container.Panel
        RecoverOriginalLoadingDirectionButton  matlab.ui.control.Button
        ThetaYEditField                 matlab.ui.control.NumericEditField
        ThetaYEditFieldLabel            matlab.ui.control.Label
        UpdateButton                    matlab.ui.control.Button
        ThetaZEditField                 matlab.ui.control.NumericEditField
        ThetaZEditFieldLabel            matlab.ui.control.Label
        ThetaXEditField                 matlab.ui.control.NumericEditField
        ThetaXEditFieldLabel            matlab.ui.control.Label
        SimulationPanel                 matlab.ui.container.Panel
        StressAnalysisonDesignButton    matlab.ui.control.Button
        RestartLinearSystemSolvingButton  matlab.ui.control.Button
        StiffnessEvaluationofVoxelbasedStructuralDesignButton  matlab.ui.control.Button
    end

    
    properties (Access = private)
        MainApp % Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            global meshHierarchy_;
            global loadingCond_;
            global loadingCondOriginal_;
            global complianceDesign_;
            global volumeFractionDesign_;
            global densityLayout_;
            app.MainApp = mainapp;
            app.DesignVolumeFractionEditField.Value = volumeFractionDesign_;
            app.DesignComplianceEditField.Value = complianceDesign_;
            if isempty(densityLayout_)
                densityLayout_ = ones(meshHierarchy_(1).numElements,1); 
            end
            app.SimulationPanel.Enable = 'off';
            loadingCondOriginal_ = loadingCond_;
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            global loadingCond_;
            global loadingCondOriginal_;
            %%Recover the original loading conditions
            loadingCond_ = loadingCondOriginal_;
            %%Reset F_
            app.ThetaXEditField.Value = 0;
            app.ThetaYEditField.Value = 0;
            app.ThetaZEditField.Value = 0; pause(0.5);
            UpdateButtonPushed(app, event);
            if isvalid(app.MainApp)
                MainWindowCtrl(app.MainApp, 1);
                app.MainApp.SimulationTasksDropDown.Value = 'None';
            end
            delete(app)            
        end

        % Button pushed function: 
        % StiffnessEvaluationofVoxelbasedStructuralDesignButton
        function StiffnessEvaluationofVoxelbasedStructuralDesignButtonPushed(app, event)
            global densityLayout_;
            global complianceDesign_;
            global volumeFractionDesign_;
            global U_; U_ = zeros(size(U_));
            
            app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel.Enable = 'off';
            app.SimulationPanel.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'off';    
            pause(1);

            GatherLSSandMPsettings(app.MainApp);
            [complianceValue, volumeFraction] = FEA_ComputeComplianceVoxel(densityLayout_);
            
            app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel.Enable = 'on';
            app.SimulationPanel.Enable = 'on';
                app.RestartLinearSystemSolvingButton.Enable = 'on';
                app.StressAnalysisonDesignButton.Enable = 'on';        
            app.ResultDisplayPanel.Enable = 'on';
                if 0==app.ThetaXEditField.Value && 0==app.ThetaYEditField.Value && 0==app.ThetaZEditField.Value
                    app.DesignComplianceEditField.Value = complianceValue;
                    app.DesignVolumeFractionEditField.Value = volumeFraction;
                    complianceDesign_ = complianceValue;
                    volumeFractionDesign_ = volumeFraction;
                    app.MainApp.DesignVolEditField.Value = volumeFractionDesign_;
                    app.MainApp.DesignComplianceEditField.Value = complianceDesign_;
                else
                    app.DesignCompliancewrtChangedLoadingDirectionsEditField.Value = complianceValue;
                end                
            ShowDeformation_Public(app.MainApp);          
        end

        % Button pushed function: RestartLinearSystemSolvingButton
        function RestartLinearSystemSolvingButtonPushed(app, event)
            global meshHierarchy_;
            global complianceDesign_;
            app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel.Enable = 'off';
            app.SimulationPanel.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'off';    
            pause(1);

            tStart = tic;
            GatherLSSandMPsettings(app.MainApp);
            Solving_CG_GMGS('printP_ON');
            ceList = TopOpti_ComputeUnitCompliance();
            complianceValue = meshHierarchy_(1).eleModulus*ceList;
            disp(['Re-start Linear System Costs: ', sprintf('%.f', toc(tStart)), 's']);
            
            app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel.Enable = 'on';
            app.SimulationPanel.Enable = 'on';
                app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Enable = 'on';
                app.RestartLinearSystemSolvingButton.Enable = 'on';
                app.StressAnalysisonDesignButton.Enable = 'on';                 
            app.ResultDisplayPanel.Enable = 'on';
                if 0==app.ThetaXEditField.Value && 0==app.ThetaYEditField.Value && 0==app.ThetaZEditField.Value
                    app.DesignComplianceEditField.Value = complianceValue;
                    complianceDesign_ = complianceValue;
                    app.MainApp.DesignComplianceEditField.Value = complianceDesign_;
                else
                    app.DesignCompliancewrtChangedLoadingDirectionsEditField.Value = complianceValue; 
                end                     
        end

        % Button pushed function: StressAnalysisonDesignButton
        function StressAnalysisonDesignButtonPushed(app, event)
            global outPath_;
            global densityLayout_;
            global meshHierarchy_;

            app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel.Enable = 'off';
            app.SimulationPanel.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'off';    
            pause(1);         
            
            disp('Stress Analysis on Design ...');            
            tStressAnalysis = tic;
            [cartesianStressFieldDesign, ~] = FEA_StressAnalysis(); 
            %%Compute per-element Von Mises stress and write the 2-channel
            %%Volume            
            vonMisesStressPerElement = FEA_ComputePerElementVonMisesStress(cartesianStressFieldDesign);
            vonMisesVolume = Common_ConvertPerEleVector2Volume(vonMisesStressPerElement);
            % niftiwrite(vonMisesVolume, strcat(outPath_, 'vonMisesStressDesignCLD.nii'));
            
            if 0==app.ThetaXEditField.Value && 0==app.ThetaYEditField.Value && 0==app.ThetaZEditField.Value && sum(densityLayout_)<meshHierarchy_(1).numElements
                IO_ExportDesignWithOneProperty_nii(vonMisesVolume, strcat(outPath_, 'ResultVolume_Design_vonMises.nii'));
                dominantDirDesign = Common_ExtractDominantDirectionsFromPrincipalStressDirections(cartesianStressFieldDesign);
                dominantDirSolid = niftiread(strcat(outPath_, 'dominantDirSolid.nii'));
                alignmentMetricVolumeByStressAlignment = Common_ComputeStressAlignmentDeviation(dominantDirSolid, dominantDirDesign);
                % niftiwrite(alignmentMetricVolumeByStressAlignment, strcat(outPath_, 'alignmentMetricVolume_byStress.nii'));            
                IO_ExportDesignWithOneProperty_nii(alignmentMetricVolumeByStressAlignment, strcat(outPath_, 'ResultVolume_Design_StressAlignment.nii'));
            else
                IO_ExportDesignWithOneProperty_nii(vonMisesVolume, strcat(outPath_, 'ResultVolume_Design_vonMises_CLD.nii'));
            end          
            disp(['Done with Stress Analysis (inc. extracting per-element Von Mises Stress Volume) after ', sprintf('%.f', toc(tStressAnalysis)), 's']);
            
            app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel.Enable = 'on';
            app.SimulationPanel.Enable = 'on';
            app.ResultDisplayPanel.Enable = 'on';
        end

        % Button pushed function: UpdateButton
        function UpdateButtonPushed(app, event)
            global meshHierarchy_;
            global loadingCond_;
            global loadingCondOriginal_;
            global F_;
            thetaX = app.ThetaXEditField.Value /180*pi;
            thetaY = app.ThetaYEditField.Value /180*pi;
            thetaZ = app.ThetaZEditField.Value /180*pi;
            
            %%Initialize Rotation Matrix
            Rx = [1 0 0; 0 cos(thetaX) -sin(thetaX); 0 sin(thetaX) cos(thetaX)];
            Ry = [cos(thetaY) 0 sin(thetaY); 0 1 0; -sin(thetaY) 0 cos(thetaY)];
            Rz = [cos(thetaZ) -sin(thetaZ) 0; sin(thetaZ) cos(thetaZ) 0; 0 0 1];
            RotMat = Rz * Ry * Rx;

            %%Rotate the Original Loading Directions
            if isempty(F_), FEA_ApplyBoundaryCondition(); end
            rotatedLoads = RotMat * loadingCondOriginal_(:,2:end)';
            loadingCond_(:,2:end) = rotatedLoads';
            F_ = sparse(meshHierarchy_(1).numNodes, 3);
	        F_(meshHierarchy_(1).nodesOnBoundary(loadingCond_(:,1)),:) = loadingCond_(:,2:end);
	        F_ = reshape(F_',meshHierarchy_(1).numDOFs,1);

            ShowProblemDescription_Public(app.MainApp);

            app.SimulationPanel.Enable = 'on';
                app.RestartLinearSystemSolvingButton.Enable = 'off';
                app.StressAnalysisonDesignButton.Enable = 'off';
        end

        % Button pushed function: RecoverOriginalLoadingDirectionButton
        function RecoverOriginalLoadingDirectionButtonPushed(app, event)
            global loadingCond_;
            global loadingCondOriginal_;
            loadingCond_ = loadingCondOriginal_;
            ShowProblemDescription_Public(app.MainApp);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 647 390];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create EvaluateDesignsStiffnesswrtChangedLoadingDirectionsPanel
            app.EvaluateDesignsStiffnesswrtChangedLoadingDirectionsPanel = uipanel(app.UIFigure);
            app.EvaluateDesignsStiffnesswrtChangedLoadingDirectionsPanel.Title = 'Evaluate Design''s Stiffness w.r.t. Changed Loading Directions';
            app.EvaluateDesignsStiffnesswrtChangedLoadingDirectionsPanel.Position = [0 8 642 383];

            % Create SimulationPanel
            app.SimulationPanel = uipanel(app.EvaluateDesignsStiffnesswrtChangedLoadingDirectionsPanel);
            app.SimulationPanel.Title = 'Simulation';
            app.SimulationPanel.Position = [1 124 640 129];

            % Create StiffnessEvaluationofVoxelbasedStructuralDesignButton
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton = uibutton(app.SimulationPanel, 'push');
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.ButtonPushedFcn = createCallbackFcn(app, @StiffnessEvaluationofVoxelbasedStructuralDesignButtonPushed, true);
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Position = [323 65 302 23];
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Text = 'Stiffness Evaluation of Voxel based Structural Design';

            % Create RestartLinearSystemSolvingButton
            app.RestartLinearSystemSolvingButton = uibutton(app.SimulationPanel, 'push');
            app.RestartLinearSystemSolvingButton.ButtonPushedFcn = createCallbackFcn(app, @RestartLinearSystemSolvingButtonPushed, true);
            app.RestartLinearSystemSolvingButton.Position = [31 65 180 23];
            app.RestartLinearSystemSolvingButton.Text = 'Re-start Linear System Solving';

            % Create StressAnalysisonDesignButton
            app.StressAnalysisonDesignButton = uibutton(app.SimulationPanel, 'push');
            app.StressAnalysisonDesignButton.ButtonPushedFcn = createCallbackFcn(app, @StressAnalysisonDesignButtonPushed, true);
            app.StressAnalysisonDesignButton.Position = [471 18 154 23];
            app.StressAnalysisonDesignButton.Text = 'Stress Analysis on Design';

            % Create RotateOriginalLoadingDirectionsviaEulersAnglesPanel
            app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel = uipanel(app.EvaluateDesignsStiffnesswrtChangedLoadingDirectionsPanel);
            app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel.Title = 'Rotate Original Loading Directions via Euler''s Angles';
            app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel.Position = [1 252 640 107];

            % Create ThetaXEditFieldLabel
            app.ThetaXEditFieldLabel = uilabel(app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel);
            app.ThetaXEditFieldLabel.HorizontalAlignment = 'right';
            app.ThetaXEditFieldLabel.Position = [37 53 48 22];
            app.ThetaXEditFieldLabel.Text = 'Theta-X';

            % Create ThetaXEditField
            app.ThetaXEditField = uieditfield(app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel, 'numeric');
            app.ThetaXEditField.Position = [100 53 100 22];

            % Create ThetaZEditFieldLabel
            app.ThetaZEditFieldLabel = uilabel(app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel);
            app.ThetaZEditFieldLabel.HorizontalAlignment = 'right';
            app.ThetaZEditFieldLabel.Position = [462 53 47 22];
            app.ThetaZEditFieldLabel.Text = 'Theta-Z';

            % Create ThetaZEditField
            app.ThetaZEditField = uieditfield(app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel, 'numeric');
            app.ThetaZEditField.Position = [524 53 100 22];

            % Create UpdateButton
            app.UpdateButton = uibutton(app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel, 'push');
            app.UpdateButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateButtonPushed, true);
            app.UpdateButton.Position = [523 12 100 23];
            app.UpdateButton.Text = 'Update';

            % Create ThetaYEditFieldLabel
            app.ThetaYEditFieldLabel = uilabel(app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel);
            app.ThetaYEditFieldLabel.HorizontalAlignment = 'right';
            app.ThetaYEditFieldLabel.Position = [264 53 48 22];
            app.ThetaYEditFieldLabel.Text = 'Theta-Y';

            % Create ThetaYEditField
            app.ThetaYEditField = uieditfield(app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel, 'numeric');
            app.ThetaYEditField.Position = [327 53 100 22];

            % Create RecoverOriginalLoadingDirectionButton
            app.RecoverOriginalLoadingDirectionButton = uibutton(app.RotateOriginalLoadingDirectionsviaEulersAnglesPanel, 'push');
            app.RecoverOriginalLoadingDirectionButton.ButtonPushedFcn = createCallbackFcn(app, @RecoverOriginalLoadingDirectionButtonPushed, true);
            app.RecoverOriginalLoadingDirectionButton.Position = [32 12 202 23];
            app.RecoverOriginalLoadingDirectionButton.Text = 'Recover Original Loading Direction';

            % Create ResultDisplayPanel
            app.ResultDisplayPanel = uipanel(app.EvaluateDesignsStiffnesswrtChangedLoadingDirectionsPanel);
            app.ResultDisplayPanel.Title = 'Result Display';
            app.ResultDisplayPanel.Position = [1 2 640 123];

            % Create DesignComplianceEditFieldLabel
            app.DesignComplianceEditFieldLabel = uilabel(app.ResultDisplayPanel);
            app.DesignComplianceEditFieldLabel.HorizontalAlignment = 'right';
            app.DesignComplianceEditFieldLabel.Position = [26 60 109 22];
            app.DesignComplianceEditFieldLabel.Text = 'Design Compliance';

            % Create DesignComplianceEditField
            app.DesignComplianceEditField = uieditfield(app.ResultDisplayPanel, 'numeric');
            app.DesignComplianceEditField.Position = [150 60 100 22];

            % Create DesignVolumeFractionEditFieldLabel
            app.DesignVolumeFractionEditFieldLabel = uilabel(app.ResultDisplayPanel);
            app.DesignVolumeFractionEditFieldLabel.HorizontalAlignment = 'right';
            app.DesignVolumeFractionEditFieldLabel.Position = [371 60 132 22];
            app.DesignVolumeFractionEditFieldLabel.Text = 'Design Volume Fraction';

            % Create DesignVolumeFractionEditField
            app.DesignVolumeFractionEditField = uieditfield(app.ResultDisplayPanel, 'numeric');
            app.DesignVolumeFractionEditField.Position = [518 60 100 22];

            % Create DesignCompliancewrtChangedLoadingDirectionsEditFieldLabel
            app.DesignCompliancewrtChangedLoadingDirectionsEditFieldLabel = uilabel(app.ResultDisplayPanel);
            app.DesignCompliancewrtChangedLoadingDirectionsEditFieldLabel.HorizontalAlignment = 'right';
            app.DesignCompliancewrtChangedLoadingDirectionsEditFieldLabel.Position = [21 21 304 22];
            app.DesignCompliancewrtChangedLoadingDirectionsEditFieldLabel.Text = 'Design Compliance (w.r.t. Changed Loading Directions)';

            % Create DesignCompliancewrtChangedLoadingDirectionsEditField
            app.DesignCompliancewrtChangedLoadingDirectionsEditField = uieditfield(app.ResultDisplayPanel, 'numeric');
            app.DesignCompliancewrtChangedLoadingDirectionsEditField.Position = [340 21 100 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Mtd_StructuralRobustnessExploration(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end