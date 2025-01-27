classdef Mtd_PSLsGuidedStructDesign < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        PSLsguidedStructuralDesignPanel  matlab.ui.container.Panel
        SettingsforPSLGenerationPanel  matlab.ui.container.Panel
        GeneratePSLsAloneButton        matlab.ui.control.Button
        MinorCheckBox                  matlab.ui.control.CheckBox
        MediumCheckBox                 matlab.ui.control.CheckBox
        MajorCheckBox                  matlab.ui.control.CheckBox
        ResultDisplayPanel             matlab.ui.container.Panel
        DesignVolumeFractionEditField  matlab.ui.control.NumericEditField
        DesignVolumeFractionEditField_2Label  matlab.ui.control.Label
        DesignComplianceEditField      matlab.ui.control.NumericEditField
        DesignComplianceEditField_2Label  matlab.ui.control.Label
        GenerationSimulationPanel      matlab.ui.container.Panel
        RestartLinearSystemSolvingButton  matlab.ui.control.Button
        StressAnalysisonDesignButton   matlab.ui.control.Button
        PSLsguidedInfillGenerationButton  matlab.ui.control.Button
        EvaluateStressAlignmentScaleButton  matlab.ui.control.Button
        StiffnessEvaluationofVoxelbasedStructuralDesignButton  matlab.ui.control.Button
        SettingsforMaterialLayoutConvertionPanel  matlab.ui.container.Panel
        FixedAreaEditField             matlab.ui.control.NumericEditField
        FixedAreaEditField_4Label      matlab.ui.control.Label
        LoadingAreaEditField           matlab.ui.control.NumericEditField
        LoadingAreaEditField_4Label    matlab.ui.control.Label
        EntireBoundaryEditField        matlab.ui.control.NumericEditField
        EntireBoundaryEditFieldLabel   matlab.ui.control.Label
        MaterialBudgetEditField        matlab.ui.control.NumericEditField
        MaterialBudgetEditFieldLabel   matlab.ui.control.Label
        TargetedTrajectoryThicknessEditField  matlab.ui.control.NumericEditField
        TargetedTrajectoryThicknessLabel  matlab.ui.control.Label
    end

    
    properties (Access = private)
        MainApp % Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            app.MainApp = mainapp;           
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Enable = 'off';
            app.RestartLinearSystemSolvingButton.Enable = 'off';
            app.StressAnalysisonDesignButton.Enable = 'off';
            app.EvaluateStressAlignmentScaleButton.Enable = 'off';
            PSLs_Preparation4TSV();
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            if isvalid(app.MainApp)
                MainWindowCtrl(app.MainApp, 1);
                app.MainApp.SimulationTasksDropDown.Value = 'None';
            end            
            delete(app)              
        end

        % Button pushed function: 
        % StiffnessEvaluationofVoxelbasedStructuralDesignButton
        function StiffnessEvaluationofVoxelbasedStructuralDesignButtonPushed(app, event)
            global complianceDesign_;
            global volumeFraction_;
            global densityLayout_;            
            global U_; U_ = zeros(size(U_));
            
            app.SettingsforPSLGenerationPanel.Enable = 'off';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'off';
            app.GenerationSimulationPanel.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'off';
                app.DesignComplianceEditField.Value = 0; 
            pause(1);            

            GatherLSSandMPsettings(app.MainApp);
            [complianceDesign_, volumeFraction_] = FEA_ComputeComplianceVoxel(densityLayout_);

            app.SettingsforPSLGenerationPanel.Enable = 'on';            
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';
                app.PSLsguidedInfillGenerationButton.Enable = 'on';
                app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Enable = 'on';
                app.RestartLinearSystemSolvingButton.Enable = 'on';
                app.StressAnalysisonDesignButton.Enable = 'on';
                app.EvaluateStressAlignmentScaleButton.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'on';
                app.DesignComplianceEditField.Value = complianceDesign_;
                app.DesignVolumeFractionEditField.Value = volumeFraction_;
            app.MainApp.DesignComplianceEditField.Value = complianceDesign_;
            app.MainApp.ShowDeformationMenu.Enable = 'on';
            ShowDeformation_Public(app.MainApp);             
        end

        % Button pushed function: EvaluateStressAlignmentScaleButton
        function EvaluateStressAlignmentScaleButtonPushed(app, event)
            global outPath_;
            
            dominantDirSolid = niftiread(strcat(outPath_, 'dominantDirSolid.nii'));
            dominantDirDesign = niftiread(strcat(outPath_, 'dominantDirDesign.nii'));
            disp('Compute Stress Aligment Scale between Solid and Design...');
            tStressAligmentAna = tic;
            alignmentMetricVolumeByStressAlignment = Common_ComputeStressAlignmentDeviation(dominantDirSolid, dominantDirDesign);            
            IO_ExportDesignWithOneProperty_nii(alignmentMetricVolumeByStressAlignment, strcat(outPath_, 'ResultVolume_Design_StressAlignment.nii'));
            disp(['Done with Stress Alignment Analysis after ', sprintf('%.1f', toc(tStressAligmentAna)), 's']);
        end

        % Button pushed function: PSLsguidedInfillGenerationButton
        function PSLsguidedInfillGenerationButtonPushed(app, event)
            global outPath_;            
            global volumeFractionDesign_;

            if ~(app.MajorCheckBox.Value || app.MediumCheckBox.Value || app.MinorCheckBox.Value)
                warning('No Principal Stress Direction Specified!'); return;
            end

            app.SettingsforPSLGenerationPanel.Enable = 'off';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'off';
            app.GenerationSimulationPanel.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'off';
                app.DesignComplianceEditField.Value = 0;
                app.DesignVolumeFractionEditField.Value = 0;                                 
            pause(1);

            psDirIndicator = zeros(1,3);
            if app.MajorCheckBox.Value, psDirIndicator(1) = 1; end
            if app.MediumCheckBox.Value, psDirIndicator(2) = 1; end
            if app.MinorCheckBox.Value, psDirIndicator(3) = 1; end
            targetDepositionRatio = app.MaterialBudgetEditField.Value;
            numLayerPSLs = app.TargetedTrajectoryThicknessEditField.Value;
            numLayerboundary = app.EntireBoundaryEditField.Value;
            numLayerLoads = app.LoadingAreaEditField.Value;
            numLayerFixation = app.FixedAreaEditField.Value;            
            if targetDepositionRatio<=0, return; end
            if numLayerPSLs<2, return; end            
            PSLs_GeneratePSLsGuidedInfillDesign(psDirIndicator, numLayerPSLs, targetDepositionRatio, numLayerboundary, numLayerLoads, numLayerFixation);
	        
            app.SettingsforPSLGenerationPanel.Enable = 'on';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';
                app.PSLsguidedInfillGenerationButton.Enable = 'on';
                app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Enable = 'on';
                app.RestartLinearSystemSolvingButton.Enable = 'off';
                app.StressAnalysisonDesignButton.Enable = 'off';
                app.EvaluateStressAlignmentScaleButton.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'on';
                app.DesignVolumeFractionEditField.Value = volumeFractionDesign_;
            app.MainApp.DesignVolEditField.Value = volumeFractionDesign_;
            app.MainApp.ShowPSLsMenu.Enable = 'on';
            app.MainApp.ShowVertexEdgeGraphMenu.Enable = 'on';
            app.MainApp.ShowDesignbyIsosurfaceNotrecommendedMenu.Enable = 'on';            
            pause(1);
            ShowPSLs_Public(app.MainApp);
            %%Output&Vis Design
            fileName = strcat(outPath_, 'ResultVolume_Design.nii');
	        IO_ExportDesignInVolume_Geo_nii(fileName);       
        end

        % Button pushed function: RestartLinearSystemSolvingButton
        function RestartLinearSystemSolvingButtonPushed(app, event)
            global meshHierarchy_;
            global complianceDesign_;

            app.SettingsforPSLGenerationPanel.Enable = 'off';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'off';
            app.GenerationSimulationPanel.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'off';
            pause(1); 

            tStart = tic;
            GatherLSSandMPsettings(app.MainApp);
            Solving_CG_GMGS('printP_ON');
            ceList = TopOpti_ComputeUnitCompliance();
            complianceDesign_ = meshHierarchy_(1).eleModulus*ceList;
            disp(['Re-start Linear System Costs: ', sprintf('%.f', toc(tStart)), 's']);

            app.SettingsforPSLGenerationPanel.Enable = 'on';            
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';
                app.PSLsguidedInfillGenerationButton.Enable = 'on';
                app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Enable = 'on';
                app.RestartLinearSystemSolvingButton.Enable = 'on';
                app.StressAnalysisonDesignButton.Enable = 'on';
                app.EvaluateStressAlignmentScaleButton.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'on';
                app.DesignComplianceEditField.Value = complianceDesign_;
            app.MainApp.DesignComplianceEditField.Value = complianceDesign_;    
            app.MainApp.ShowDeformationMenu.Enable = 'on';
            ShowDeformation_Public(app.MainApp);    
        end

        % Button pushed function: StressAnalysisonDesignButton
        function StressAnalysisonDesignButtonPushed(app, event)
            global outPath_;            
            app.SettingsforPSLGenerationPanel.Enable = 'off';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'off';
            app.GenerationSimulationPanel.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'off';
            pause(1); 
           
            disp('Stress Analysis on Design ...');            
            tStressAnalysis = tic;
            [cartesianStressFieldDesign, ~] = FEA_StressAnalysis();  
            vonMisesStressPerElement = FEA_ComputePerElementVonMisesStress(cartesianStressFieldDesign);
            dominantDirDesign = Common_ExtractDominantDirectionsFromPrincipalStressDirections(cartesianStressFieldDesign);
            if ~isempty(dominantDirDesign)
                niftiwrite(dominantDirDesign, strcat(outPath_, 'dominantDirDesign.nii'));
            end
            vonMisesVolume = Common_ConvertPerEleVector2Volume(vonMisesStressPerElement);
            IO_ExportDesignWithOneProperty_nii(vonMisesVolume, strcat(outPath_, 'ResultVolume_Design_vonMises.nii'));              
            disp(['Done with Stress Analysis (inc. extracting dominant stress directions) after ', sprintf('%.f', toc(tStressAnalysis)), 's']);

            app.SettingsforPSLGenerationPanel.Enable = 'on';            
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';
                app.PSLsguidedInfillGenerationButton.Enable = 'on';
                app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Enable = 'on';
                app.RestartLinearSystemSolvingButton.Enable = 'on';
                app.StressAnalysisonDesignButton.Enable = 'on';
                app.EvaluateStressAlignmentScaleButton.Enable = 'on';
            app.ResultDisplayPanel.Enable = 'on';            
        end

        % Button pushed function: GeneratePSLsAloneButton
        function GeneratePSLsAloneButtonPushed(app, event)
            if ~(app.MajorCheckBox.Value || app.MediumCheckBox.Value || app.MinorCheckBox.Value)
                warning('No Principal Stress Direction Specified!'); return;
            end
            
            app.SettingsforPSLGenerationPanel.Enable = 'off';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'off';
            app.GenerationSimulationPanel.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'off';                            
            pause(1);

            psDirIndicator = zeros(1,3);
            if app.MajorCheckBox.Value, psDirIndicator(1) = 1; end
            if app.MediumCheckBox.Value, psDirIndicator(2) = 1; end
            if app.MinorCheckBox.Value, psDirIndicator(3) = 1; end

            lineDensCtrl = 10;
            PSLs_GeneratePSLsBy3DTSV(lineDensCtrl, psDirIndicator);

            app.SettingsforPSLGenerationPanel.Enable = 'on';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';
            app.ResultDisplayPanel.Enable = 'on';
            app.MainApp.ShowPSLsMenu.Enable = 'on';          
            pause(1);
            ShowPSLs_Public(app.MainApp);            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 654 586];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create PSLsguidedStructuralDesignPanel
            app.PSLsguidedStructuralDesignPanel = uipanel(app.UIFigure);
            app.PSLsguidedStructuralDesignPanel.Title = 'PSLs-guided Structural Design';
            app.PSLsguidedStructuralDesignPanel.Position = [0 19 641 568];

            % Create SettingsforMaterialLayoutConvertionPanel
            app.SettingsforMaterialLayoutConvertionPanel = uipanel(app.PSLsguidedStructuralDesignPanel);
            app.SettingsforMaterialLayoutConvertionPanel.Title = 'Settings for Material Layout Convertion';
            app.SettingsforMaterialLayoutConvertionPanel.Position = [0 281 640 169];

            % Create TargetedTrajectoryThicknessLabel
            app.TargetedTrajectoryThicknessLabel = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.TargetedTrajectoryThicknessLabel.HorizontalAlignment = 'right';
            app.TargetedTrajectoryThicknessLabel.Position = [379 57 165 22];
            app.TargetedTrajectoryThicknessLabel.Text = 'Targeted Trajectory Thickness';

            % Create TargetedTrajectoryThicknessEditField
            app.TargetedTrajectoryThicknessEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.TargetedTrajectoryThicknessEditField.ValueDisplayFormat = '%.0f';
            app.TargetedTrajectoryThicknessEditField.Position = [559 57 68 22];
            app.TargetedTrajectoryThicknessEditField.Value = 3;

            % Create MaterialBudgetEditFieldLabel
            app.MaterialBudgetEditFieldLabel = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.MaterialBudgetEditFieldLabel.HorizontalAlignment = 'right';
            app.MaterialBudgetEditFieldLabel.Position = [455 17 89 22];
            app.MaterialBudgetEditFieldLabel.Text = 'Material Budget';

            % Create MaterialBudgetEditField
            app.MaterialBudgetEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.MaterialBudgetEditField.Position = [559 17 68 22];
            app.MaterialBudgetEditField.Value = 0.4;

            % Create EntireBoundaryEditFieldLabel
            app.EntireBoundaryEditFieldLabel = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.EntireBoundaryEditFieldLabel.HorizontalAlignment = 'right';
            app.EntireBoundaryEditFieldLabel.Position = [453 100 91 22];
            app.EntireBoundaryEditFieldLabel.Text = 'Entire Boundary';

            % Create EntireBoundaryEditField
            app.EntireBoundaryEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.EntireBoundaryEditField.ValueDisplayFormat = '%.0f';
            app.EntireBoundaryEditField.Position = [559 100 68 22];
            app.EntireBoundaryEditField.Value = 2;

            % Create LoadingAreaEditField_4Label
            app.LoadingAreaEditField_4Label = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.LoadingAreaEditField_4Label.HorizontalAlignment = 'right';
            app.LoadingAreaEditField_4Label.Position = [39 100 76 22];
            app.LoadingAreaEditField_4Label.Text = 'Loading Area';

            % Create LoadingAreaEditField
            app.LoadingAreaEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.LoadingAreaEditField.ValueDisplayFormat = '%.0f';
            app.LoadingAreaEditField.Position = [130 100 68 22];

            % Create FixedAreaEditField_4Label
            app.FixedAreaEditField_4Label = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.FixedAreaEditField_4Label.HorizontalAlignment = 'right';
            app.FixedAreaEditField_4Label.Position = [267 100 62 22];
            app.FixedAreaEditField_4Label.Text = 'Fixed Area';

            % Create FixedAreaEditField
            app.FixedAreaEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.FixedAreaEditField.ValueDisplayFormat = '%.0f';
            app.FixedAreaEditField.Position = [344 100 68 22];

            % Create GenerationSimulationPanel
            app.GenerationSimulationPanel = uipanel(app.PSLsguidedStructuralDesignPanel);
            app.GenerationSimulationPanel.Title = 'Generation & Simulation';
            app.GenerationSimulationPanel.Position = [0 82 640 200];

            % Create StiffnessEvaluationofVoxelbasedStructuralDesignButton
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.ButtonPushedFcn = createCallbackFcn(app, @StiffnessEvaluationofVoxelbasedStructuralDesignButtonPushed, true);
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Position = [325 101 302 23];
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Text = 'Stiffness Evaluation of Voxel based Structural Design';

            % Create EvaluateStressAlignmentScaleButton
            app.EvaluateStressAlignmentScaleButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.EvaluateStressAlignmentScaleButton.ButtonPushedFcn = createCallbackFcn(app, @EvaluateStressAlignmentScaleButtonPushed, true);
            app.EvaluateStressAlignmentScaleButton.Position = [439 18 188 23];
            app.EvaluateStressAlignmentScaleButton.Text = 'Evaluate Stress Alignment Scale';

            % Create PSLsguidedInfillGenerationButton
            app.PSLsguidedInfillGenerationButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.PSLsguidedInfillGenerationButton.ButtonPushedFcn = createCallbackFcn(app, @PSLsguidedInfillGenerationButtonPushed, true);
            app.PSLsguidedInfillGenerationButton.BackgroundColor = [0.9608 0.9608 0.9608];
            app.PSLsguidedInfillGenerationButton.FontWeight = 'bold';
            app.PSLsguidedInfillGenerationButton.Position = [444 142 183 23];
            app.PSLsguidedInfillGenerationButton.Text = 'PSLs-guided Infill Generation';

            % Create StressAnalysisonDesignButton
            app.StressAnalysisonDesignButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.StressAnalysisonDesignButton.ButtonPushedFcn = createCallbackFcn(app, @StressAnalysisonDesignButtonPushed, true);
            app.StressAnalysisonDesignButton.Position = [473 59 154 23];
            app.StressAnalysisonDesignButton.Text = 'Stress Analysis on Design';

            % Create RestartLinearSystemSolvingButton
            app.RestartLinearSystemSolvingButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.RestartLinearSystemSolvingButton.ButtonPushedFcn = createCallbackFcn(app, @RestartLinearSystemSolvingButtonPushed, true);
            app.RestartLinearSystemSolvingButton.Position = [117 101 180 23];
            app.RestartLinearSystemSolvingButton.Text = 'Re-start Linear System Solving';

            % Create ResultDisplayPanel
            app.ResultDisplayPanel = uipanel(app.PSLsguidedStructuralDesignPanel);
            app.ResultDisplayPanel.Title = 'Result Display';
            app.ResultDisplayPanel.Position = [0 0 640 83];

            % Create DesignComplianceEditField_2Label
            app.DesignComplianceEditField_2Label = uilabel(app.ResultDisplayPanel);
            app.DesignComplianceEditField_2Label.HorizontalAlignment = 'right';
            app.DesignComplianceEditField_2Label.Position = [85 21 109 22];
            app.DesignComplianceEditField_2Label.Text = 'Design Compliance';

            % Create DesignComplianceEditField
            app.DesignComplianceEditField = uieditfield(app.ResultDisplayPanel, 'numeric');
            app.DesignComplianceEditField.Position = [209 21 100 22];

            % Create DesignVolumeFractionEditField_2Label
            app.DesignVolumeFractionEditField_2Label = uilabel(app.ResultDisplayPanel);
            app.DesignVolumeFractionEditField_2Label.HorizontalAlignment = 'right';
            app.DesignVolumeFractionEditField_2Label.Position = [380 21 132 22];
            app.DesignVolumeFractionEditField_2Label.Text = 'Design Volume Fraction';

            % Create DesignVolumeFractionEditField
            app.DesignVolumeFractionEditField = uieditfield(app.ResultDisplayPanel, 'numeric');
            app.DesignVolumeFractionEditField.Position = [527 21 100 22];

            % Create SettingsforPSLGenerationPanel
            app.SettingsforPSLGenerationPanel = uipanel(app.PSLsguidedStructuralDesignPanel);
            app.SettingsforPSLGenerationPanel.Title = 'Settings for PSL Generation';
            app.SettingsforPSLGenerationPanel.Position = [1 450 640 96];

            % Create MajorCheckBox
            app.MajorCheckBox = uicheckbox(app.SettingsforPSLGenerationPanel);
            app.MajorCheckBox.Text = 'Major';
            app.MajorCheckBox.Position = [15 27 52 22];
            app.MajorCheckBox.Value = true;

            % Create MediumCheckBox
            app.MediumCheckBox = uicheckbox(app.SettingsforPSLGenerationPanel);
            app.MediumCheckBox.Text = 'Medium';
            app.MediumCheckBox.Position = [168 27 65 22];

            % Create MinorCheckBox
            app.MinorCheckBox = uicheckbox(app.SettingsforPSLGenerationPanel);
            app.MinorCheckBox.Text = 'Minor';
            app.MinorCheckBox.Position = [345 27 52 22];
            app.MinorCheckBox.Value = true;

            % Create GeneratePSLsAloneButton
            app.GeneratePSLsAloneButton = uibutton(app.SettingsforPSLGenerationPanel, 'push');
            app.GeneratePSLsAloneButton.ButtonPushedFcn = createCallbackFcn(app, @GeneratePSLsAloneButtonPushed, true);
            app.GeneratePSLsAloneButton.Position = [497 27 130 23];
            app.GeneratePSLsAloneButton.Text = 'Generate PSLs Alone';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Mtd_PSLsGuidedStructDesign(varargin)

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