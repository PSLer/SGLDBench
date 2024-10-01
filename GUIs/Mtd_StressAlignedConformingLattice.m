classdef Mtd_StressAlignedConformingLattice < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        StressalignedGeometricParametrizationbasedDesignPanel  matlab.ui.container.Panel
        InputDataPreparationTetmeshStressFieldPanel  matlab.ui.container.Panel
        TetrahedraApproxiEditField      matlab.ui.control.NumericEditField
        TetrahedraApproxiEditFieldLabel  matlab.ui.control.Label
        DataPreparationButton           matlab.ui.control.Button
        RefTextArea                     matlab.ui.control.TextArea
        Label                           matlab.ui.control.Label
        ResultDisplayPanel              matlab.ui.container.Panel
        DesignVolumeFractionEditField   matlab.ui.control.NumericEditField
        DesignVolumeFractionEditFieldLabel  matlab.ui.control.Label
        DesignComplianceEditField       matlab.ui.control.NumericEditField
        DesignComplianceEditFieldLabel  matlab.ui.control.Label
        GenerationSimulationPanel       matlab.ui.container.Panel
        LatticeAspectRatioEditField     matlab.ui.control.NumericEditField
        LatticeAspectRatioEditFieldLabel  matlab.ui.control.Label
        EvaluateStressAlignmentScaleButton  matlab.ui.control.Button
        StressAnalysisonDesignButton    matlab.ui.control.Button
        RestartLinearSystemSolvingButton  matlab.ui.control.Button
        StressalignedConformingLatticeInfillDesignButton  matlab.ui.control.Button
        StiffnessEvaluationofVoxelbasedStructuralDesignButton  matlab.ui.control.Button
        SettingsforMaterialLayoutConvertionPanel  matlab.ui.container.Panel
        MaterialBudgetEditField         matlab.ui.control.NumericEditField
        MaterialBudgetEditFieldLabel    matlab.ui.control.Label
        FixedAreaEditField              matlab.ui.control.NumericEditField
        FixedAreaEditFieldLabel         matlab.ui.control.Label
        LoadingAreaEditField            matlab.ui.control.NumericEditField
        LoadingAreaEditFieldLabel       matlab.ui.control.Label
        EntireBoundaryEditField         matlab.ui.control.NumericEditField
        EntireBoundaryEditFieldLabel    matlab.ui.control.Label
        EdgeThicknessEditField          matlab.ui.control.NumericEditField
        EdgeThicknessEditFieldLabel     matlab.ui.control.Label
    end

    
    properties (Access = private)
        MainApp % Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            addpath('./externalModules/Gao2017/');
            app.MainApp = mainapp;
            app.GenerationSimulationPanel.Enable = 'off';
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            if isvalid(app.MainApp)
                MainWindowCtrl(app.MainApp, 1);
                app.MainApp.SimulationTasksDropDown.Value = 'None';
            end
            delete(app) 
        end

        % Button pushed function: DataPreparationButton
        function DataPreparationButtonPushed(app, event)
            global cartesianStressField_;
            app.DataPreparationButton.Enable = 'off';
            app.GenerationSimulationPanel.Enable = 'off';
            pause(1);

            if isempty(cartesianStressField_)
                warning('No Stress Field is Available!');
                app.DataPreparationButton.Enable = 'on';
                return;
            end
            
            disp('Input Data is Preparing...');
            numTetrahedraInGatewayMesh = app.TetrahedraApproxiEditField.Value;
            numTetrahedraInGatewayMesh = max(numTetrahedraInGatewayMesh, 50000);
            SAGS_GenerateDelaunayTetMeshFromInputSurfaceMesh(round(numTetrahedraInGatewayMesh/2));            
            SAGS_InterpolatingStressFieldOnTetMesh();

            app.DataPreparationButton.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';
                app.StressalignedConformingLatticeInfillDesignButton.Enable = 'on';
                app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Enable = 'off';
                app.RestartLinearSystemSolvingButton.Enable = 'off';
                app.StressAnalysisonDesignButton.Enable = 'off';
                app.EvaluateStressAlignmentScaleButton.Enable = 'off';                    
            disp('Input Data is Ready!');
        end

        % Button pushed function: 
        % StressalignedConformingLatticeInfillDesignButton
        function StressalignedConformingLatticeInfillDesignButtonPushed(app, event)
            global outPath_;
            global volumeFractionDesign_;

            app.InputDataPreparationTetmeshStressFieldPanel.Enable = 'off';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'off';
            app.GenerationSimulationPanel.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'off';
                app.DesignComplianceEditField.Value = 0;
                app.DesignVolumeFractionEditField.Value = 0;            
            pause(1);

            targetDepositionRatio = app.MaterialBudgetEditField.Value;
            edgeWidth = app.EdgeThicknessEditField.Value; edgeWidth = max(edgeWidth,1);
            numLayerboundary = app.EntireBoundaryEditField.Value; numLayerboundary = max(numLayerboundary,0);
            numLayerLoads = app.LoadingAreaEditField.Value; numLayerLoads = max(numLayerLoads,0);
            numLayerFixation = app.FixedAreaEditField.Value; numLayerFixation = max(numLayerFixation,0);
            aspectRatio = app.LatticeAspectRatioEditField.Value;
            if aspectRatio<0, return; end
            aspectRatio = min(aspectRatio,1);
            aspectRatio = max(aspectRatio,0.1);
            app.LatticeAspectRatioEditField.Value = aspectRatio;
            if targetDepositionRatio<=0, return; end
            
            SAGS_StressAlignedConformingLatticeGeneration(edgeWidth, targetDepositionRatio, numLayerboundary, ...
                numLayerLoads, numLayerFixation, aspectRatio);
            
            app.InputDataPreparationTetmeshStressFieldPanel.Enable = 'on';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';
                app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Enable = 'on';
            app.ResultDisplayPanel.Enable = 'on';
                app.DesignVolumeFractionEditField.Value = volumeFractionDesign_;    
            app.MainApp.DesignVolEditField.Value = volumeFractionDesign_;
            app.MainApp.ShowVertexEdgeGraphMenu.Enable = 'on';
            app.MainApp.ShowDesignbyDensityFieldNotrecommendedMenu.Enable = 'on';
            ShowEdgeVertexGraph_Public(app.MainApp);            

            %%Output&Vis Design
            fileName = strcat(outPath_, 'DesignVolume.nii');
	        IO_ExportDesignInVolume_Geo_nii(fileName);
            %system('"./src/quokka.exe" ./out/DesignVolume.nii');    
        end

        % Button pushed function: 
        % StiffnessEvaluationofVoxelbasedStructuralDesignButton
        function StiffnessEvaluationofVoxelbasedStructuralDesignButtonPushed(app, event)
            global complianceDesign_;
            global volumeFraction_;
            global densityLayout_;            
            global U_; U_ = zeros(size(U_));
            
            app.InputDataPreparationTetmeshStressFieldPanel.Enable = 'off';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'off';
            app.GenerationSimulationPanel.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'off';    
            pause(1);

            GatherLSSandMPsettings(app.MainApp);
            [complianceDesign_, volumeFraction_] = FEA_ComputeComplianceVoxel(densityLayout_);
            
            app.InputDataPreparationTetmeshStressFieldPanel.Enable = 'on';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';
                app.RestartLinearSystemSolvingButton.Enable = 'on';
                app.StressAnalysisonDesignButton.Enable = 'on';                
            app.ResultDisplayPanel.Enable = 'on'; 
                app.DesignComplianceEditField.Value = complianceDesign_;
            app.MainApp.DesignComplianceEditField.Value = complianceDesign_;
            % app.MainApp.SimulationTasksDropDown.Enable = 'on';
            % app.MainApp.ShowDeformationMenu.Enable = 'on';
            ShowDeformation_Public(app.MainApp);            
        end

        % Button pushed function: RestartLinearSystemSolvingButton
        function RestartLinearSystemSolvingButtonPushed(app, event)
            global meshHierarchy_;
            global complianceDesign_;

            app.InputDataPreparationTetmeshStressFieldPanel.Enable = 'off';
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
            
            app.InputDataPreparationTetmeshStressFieldPanel.Enable = 'on';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';              
            app.ResultDisplayPanel.Enable = 'on'; 
                app.DesignComplianceEditField.Value = complianceDesign_;
            app.MainApp.DesignComplianceEditField.Value = complianceDesign_;    
        end

        % Button pushed function: StressAnalysisonDesignButton
        function StressAnalysisonDesignButtonPushed(app, event)
            global outPath_;            

            app.InputDataPreparationTetmeshStressFieldPanel.Enable = 'off';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'off';
            app.GenerationSimulationPanel.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'off';    
            pause(1);         
            
            disp('Stress Analysis on Design ...');            
            tStressAnalysis = tic;
            [cartesianStressFieldDesign, ~] = FEA_StressAnalysis();  
            dominantDirDesign = Common_ExtractDominantDirectionsFromPrincipalStressDirections(cartesianStressFieldDesign);             
            if ~isempty(dominantDirDesign)
                niftiwrite(dominantDirDesign, strcat(outPath_, 'dominantDirDesign.nii'));
            end
            disp(['Done with Stress Analysis (inc. extracting dominant stress directions) after ', sprintf('%.f', toc(tStressAnalysis)), 's']);
            
            app.InputDataPreparationTetmeshStressFieldPanel.Enable = 'on';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';
                app.EvaluateStressAlignmentScaleButton.Enable = 'on';
            app.ResultDisplayPanel.Enable = 'on';             
        end

        % Button pushed function: EvaluateStressAlignmentScaleButton
        function EvaluateStressAlignmentScaleButtonPushed(app, event)
            global outPath_;
            app.EvaluateStressAlignmentScaleButton.Enable = 'off'; 
            pause(1);

            dominantDirSolid = niftiread(strcat(outPath_, 'dominantDirSolid.nii'));
            dominantDirDesign = niftiread(strcat(outPath_, 'dominantDirDesign.nii'));
            disp('Compute Stress Aligment Scale between Solid and Design...');
            tStressAligmentAna = tic;
            alignmentMetricVolumeByStressAlignment = Common_ComputeStressAlignmentDeviation(dominantDirSolid, dominantDirDesign);
            niftiwrite(alignmentMetricVolumeByStressAlignment, strcat(outPath_, 'alignmentMetricVolume_byStress.nii'));            
            % alignmentMetricVolumeByEdgeAlignment = Common_ComputeEdgeAlignmentDeviation(dominantDirDesign);
            % niftiwrite(alignmentMetricVolumeByEdgeAlignment, strcat(outPath_, 'alignmentMetricVolume_byEdge.nii'));
            disp(['Done with Stress Alignment Analysis after ', sprintf('%.1f', toc(tStressAligmentAna)), 's']);
            app.EvaluateStressAlignmentScaleButton.Enable = 'on'; 
            %system('"./src/quokka.exe" ./out/alignmentMetricVolume_byStress.nii');                
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 654 651];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create StressalignedGeometricParametrizationbasedDesignPanel
            app.StressalignedGeometricParametrizationbasedDesignPanel = uipanel(app.UIFigure);
            app.StressalignedGeometricParametrizationbasedDesignPanel.Title = 'Conforming Lattice Structure';
            app.StressalignedGeometricParametrizationbasedDesignPanel.Position = [0 13 640 639];

            % Create SettingsforMaterialLayoutConvertionPanel
            app.SettingsforMaterialLayoutConvertionPanel = uipanel(app.StressalignedGeometricParametrizationbasedDesignPanel);
            app.SettingsforMaterialLayoutConvertionPanel.Title = 'Settings for Material Layout Convertion';
            app.SettingsforMaterialLayoutConvertionPanel.Position = [0 404 640 142];

            % Create EdgeThicknessEditFieldLabel
            app.EdgeThicknessEditFieldLabel = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.EdgeThicknessEditFieldLabel.HorizontalAlignment = 'right';
            app.EdgeThicknessEditFieldLabel.Position = [453 49 90 22];
            app.EdgeThicknessEditFieldLabel.Text = 'Edge Thickness';

            % Create EdgeThicknessEditField
            app.EdgeThicknessEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.EdgeThicknessEditField.ValueDisplayFormat = '%.0f';
            app.EdgeThicknessEditField.Position = [558 49 68 22];
            app.EdgeThicknessEditField.Value = 3;

            % Create EntireBoundaryEditFieldLabel
            app.EntireBoundaryEditFieldLabel = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.EntireBoundaryEditFieldLabel.HorizontalAlignment = 'right';
            app.EntireBoundaryEditFieldLabel.Position = [451 87 91 22];
            app.EntireBoundaryEditFieldLabel.Text = 'Entire Boundary';

            % Create EntireBoundaryEditField
            app.EntireBoundaryEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.EntireBoundaryEditField.ValueDisplayFormat = '%.0f';
            app.EntireBoundaryEditField.Position = [557 87 68 22];
            app.EntireBoundaryEditField.Value = 2;

            % Create LoadingAreaEditFieldLabel
            app.LoadingAreaEditFieldLabel = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.LoadingAreaEditFieldLabel.HorizontalAlignment = 'right';
            app.LoadingAreaEditFieldLabel.Position = [65 87 76 22];
            app.LoadingAreaEditFieldLabel.Text = 'Loading Area';

            % Create LoadingAreaEditField
            app.LoadingAreaEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.LoadingAreaEditField.ValueDisplayFormat = '%.0f';
            app.LoadingAreaEditField.Position = [156 87 68 22];

            % Create FixedAreaEditFieldLabel
            app.FixedAreaEditFieldLabel = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.FixedAreaEditFieldLabel.HorizontalAlignment = 'right';
            app.FixedAreaEditFieldLabel.Position = [272 87 62 22];
            app.FixedAreaEditFieldLabel.Text = 'Fixed Area';

            % Create FixedAreaEditField
            app.FixedAreaEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.FixedAreaEditField.ValueDisplayFormat = '%.0f';
            app.FixedAreaEditField.Position = [349 87 68 22];

            % Create MaterialBudgetEditFieldLabel
            app.MaterialBudgetEditFieldLabel = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.MaterialBudgetEditFieldLabel.HorizontalAlignment = 'right';
            app.MaterialBudgetEditFieldLabel.Position = [454 11 89 22];
            app.MaterialBudgetEditFieldLabel.Text = 'Material Budget';

            % Create MaterialBudgetEditField
            app.MaterialBudgetEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.MaterialBudgetEditField.Position = [558 11 68 22];
            app.MaterialBudgetEditField.Value = 0.4;

            % Create GenerationSimulationPanel
            app.GenerationSimulationPanel = uipanel(app.StressalignedGeometricParametrizationbasedDesignPanel);
            app.GenerationSimulationPanel.Title = 'Generation & Simulation';
            app.GenerationSimulationPanel.Position = [0 191 640 213];

            % Create StiffnessEvaluationofVoxelbasedStructuralDesignButton
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.ButtonPushedFcn = createCallbackFcn(app, @StiffnessEvaluationofVoxelbasedStructuralDesignButtonPushed, true);
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Position = [325 108 302 23];
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Text = 'Stiffness Evaluation of Voxel based Structural Design';

            % Create StressalignedConformingLatticeInfillDesignButton
            app.StressalignedConformingLatticeInfillDesignButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.StressalignedConformingLatticeInfillDesignButton.ButtonPushedFcn = createCallbackFcn(app, @StressalignedConformingLatticeInfillDesignButtonPushed, true);
            app.StressalignedConformingLatticeInfillDesignButton.BackgroundColor = [0.9608 0.9608 0.9608];
            app.StressalignedConformingLatticeInfillDesignButton.Position = [363 152 264 23];
            app.StressalignedConformingLatticeInfillDesignButton.Text = 'Stress-aligned Conforming Lattice Infill Design';

            % Create RestartLinearSystemSolvingButton
            app.RestartLinearSystemSolvingButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.RestartLinearSystemSolvingButton.ButtonPushedFcn = createCallbackFcn(app, @RestartLinearSystemSolvingButtonPushed, true);
            app.RestartLinearSystemSolvingButton.Position = [113 108 180 23];
            app.RestartLinearSystemSolvingButton.Text = 'Re-start Linear System Solving';

            % Create StressAnalysisonDesignButton
            app.StressAnalysisonDesignButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.StressAnalysisonDesignButton.ButtonPushedFcn = createCallbackFcn(app, @StressAnalysisonDesignButtonPushed, true);
            app.StressAnalysisonDesignButton.Position = [473 62 154 23];
            app.StressAnalysisonDesignButton.Text = 'Stress Analysis on Design';

            % Create EvaluateStressAlignmentScaleButton
            app.EvaluateStressAlignmentScaleButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.EvaluateStressAlignmentScaleButton.ButtonPushedFcn = createCallbackFcn(app, @EvaluateStressAlignmentScaleButtonPushed, true);
            app.EvaluateStressAlignmentScaleButton.Position = [438 17 188 23];
            app.EvaluateStressAlignmentScaleButton.Text = 'Evaluate Stress Alignment Scale';

            % Create LatticeAspectRatioEditFieldLabel
            app.LatticeAspectRatioEditFieldLabel = uilabel(app.GenerationSimulationPanel);
            app.LatticeAspectRatioEditFieldLabel.HorizontalAlignment = 'right';
            app.LatticeAspectRatioEditFieldLabel.Position = [81 152 111 22];
            app.LatticeAspectRatioEditFieldLabel.Text = 'Lattice Aspect Ratio';

            % Create LatticeAspectRatioEditField
            app.LatticeAspectRatioEditField = uieditfield(app.GenerationSimulationPanel, 'numeric');
            app.LatticeAspectRatioEditField.Position = [207 152 100 22];
            app.LatticeAspectRatioEditField.Value = 1;

            % Create ResultDisplayPanel
            app.ResultDisplayPanel = uipanel(app.StressalignedGeometricParametrizationbasedDesignPanel);
            app.ResultDisplayPanel.Title = 'Result Display';
            app.ResultDisplayPanel.Position = [0 88 640 104];

            % Create DesignComplianceEditFieldLabel
            app.DesignComplianceEditFieldLabel = uilabel(app.ResultDisplayPanel);
            app.DesignComplianceEditFieldLabel.HorizontalAlignment = 'right';
            app.DesignComplianceEditFieldLabel.Position = [26 30 109 22];
            app.DesignComplianceEditFieldLabel.Text = 'Design Compliance';

            % Create DesignComplianceEditField
            app.DesignComplianceEditField = uieditfield(app.ResultDisplayPanel, 'numeric');
            app.DesignComplianceEditField.Position = [150 30 100 22];

            % Create DesignVolumeFractionEditFieldLabel
            app.DesignVolumeFractionEditFieldLabel = uilabel(app.ResultDisplayPanel);
            app.DesignVolumeFractionEditFieldLabel.HorizontalAlignment = 'right';
            app.DesignVolumeFractionEditFieldLabel.Position = [379 30 132 22];
            app.DesignVolumeFractionEditFieldLabel.Text = 'Design Volume Fraction';

            % Create DesignVolumeFractionEditField
            app.DesignVolumeFractionEditField = uieditfield(app.ResultDisplayPanel, 'numeric');
            app.DesignVolumeFractionEditField.Position = [526 30 100 22];

            % Create Label
            app.Label = uilabel(app.StressalignedGeometricParametrizationbasedDesignPanel);
            app.Label.HorizontalAlignment = 'right';
            app.Label.Position = [5 56 27 22];
            app.Label.Text = 'Ref.';

            % Create RefTextArea
            app.RefTextArea = uitextarea(app.StressalignedGeometricParametrizationbasedDesignPanel);
            app.RefTextArea.Position = [46 10 594 70];
            app.RefTextArea.Value = {'[1] Gao, Xifeng, et al. "Robust hex-dominant mesh generation using field-guided polyhedral agglomeration." ACM Transactions on Graphics (TOG) 36.4 (2017): 1-13.'; '[2] Wu, Jun, Weiming Wang, and Xifeng Gao. "Design and optimization of conforming lattice structures." IEEE transactions on visualization and computer graphics 27.1 (2019): 43-56.'};

            % Create InputDataPreparationTetmeshStressFieldPanel
            app.InputDataPreparationTetmeshStressFieldPanel = uipanel(app.StressalignedGeometricParametrizationbasedDesignPanel);
            app.InputDataPreparationTetmeshStressFieldPanel.Title = 'Input Data Preparation (Tet-mesh & Stress Field)';
            app.InputDataPreparationTetmeshStressFieldPanel.Position = [0 546 640 71];

            % Create DataPreparationButton
            app.DataPreparationButton = uibutton(app.InputDataPreparationTetmeshStressFieldPanel, 'push');
            app.DataPreparationButton.ButtonPushedFcn = createCallbackFcn(app, @DataPreparationButtonPushed, true);
            app.DataPreparationButton.Position = [521 13 106 23];
            app.DataPreparationButton.Text = 'Data Preparation';

            % Create TetrahedraApproxiEditFieldLabel
            app.TetrahedraApproxiEditFieldLabel = uilabel(app.InputDataPreparationTetmeshStressFieldPanel);
            app.TetrahedraApproxiEditFieldLabel.HorizontalAlignment = 'right';
            app.TetrahedraApproxiEditFieldLabel.Position = [104 13 116 22];
            app.TetrahedraApproxiEditFieldLabel.Text = '#Tetrahedra Approxi.';

            % Create TetrahedraApproxiEditField
            app.TetrahedraApproxiEditField = uieditfield(app.InputDataPreparationTetmeshStressFieldPanel, 'numeric');
            app.TetrahedraApproxiEditField.ValueDisplayFormat = '%.0f';
            app.TetrahedraApproxiEditField.Position = [235 13 68 22];
            app.TetrahedraApproxiEditField.Value = 500000;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Mtd_StressAlignedConformingLattice(varargin)

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