classdef Mtd_MeshGraphBasedStructDesign < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        OpenMenu                        matlab.ui.container.Menu
        SolidHexTetMeshmeshvtkmshMenu   matlab.ui.container.Menu
        GraphobjMenu                    matlab.ui.container.Menu
        EvaluateExternalMeshGraphbasedStructuralDesignPanel  matlab.ui.container.Panel
        ResultDisplayPanel              matlab.ui.container.Panel
        DesignVolumeFractionEditField   matlab.ui.control.NumericEditField
        DesignVolumeFractionEditFieldLabel  matlab.ui.control.Label
        DesignComplianceEditField       matlab.ui.control.NumericEditField
        DesignComplianceEditFieldLabel  matlab.ui.control.Label
        PreProcessPanel                 matlab.ui.container.Panel
        DataAlignmentButton             matlab.ui.control.Button
        GenerationSimulationPanel       matlab.ui.container.Panel
        GenerateVoxelbasedStructuralDesignButton  matlab.ui.control.Button
        StressAnalysisonDesignButton    matlab.ui.control.Button
        RestartLinearSystemSolvingButton  matlab.ui.control.Button
        EvaluateStressAlignmentScaleButton  matlab.ui.control.Button
        StiffnessEvaluationofVoxelbasedStructuralDesignButton  matlab.ui.control.Button
        SettingsforMaterialLayoutConvertionPanel  matlab.ui.container.Panel
        FixedAreaEditField              matlab.ui.control.NumericEditField
        FixedAreaEditField_4Label       matlab.ui.control.Label
        LoadingAreaEditField            matlab.ui.control.NumericEditField
        LoadingAreaEditField_4Label     matlab.ui.control.Label
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
            app.MainApp = mainapp;            
            app.PreProcessPanel.Enable = 'off';
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

        % Button pushed function: DataAlignmentButton
        function DataAlignmentButtonPushed(app, event)
            global axHandle_;
            global meshHierarchy_;
            global frameStruct4Voxelization_;

            app.PreProcessPanel.Enable = 'off';
            app.GenerationSimulationPanel.Enable = 'off';
            pause(1);

            MGD_DataPreprocess();
            
            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            [az, el] = view(axHandle_);
            cla(axHandle_); colorbar(axHandle_, 'off');
            hdSilhouette = Vis_DrawMesh3D(axHandle_, meshHierarchy_(1).boundaryNodeCoords, meshHierarchy_(1).boundaryEleFaces, 0);
            set(hdSilhouette, 'FaceColor', [65 174 118]/255, 'FaceAlpha', 0.15);
            hold(axHandle_, 'on');
            Vis_DrawGraph3D(axHandle_, frameStruct4Voxelization_.nodeCoords, frameStruct4Voxelization_.eNodMat, 500);
            view(axHandle_, az, el);
            Vis_UserLighting(axHandle_);

            app.PreProcessPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';
                app.GenerateVoxelbasedStructuralDesignButton.Enable = 'on';
                app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Enable = 'off';
                app.RestartLinearSystemSolvingButton.Enable = 'off';
                app.StressAnalysisonDesignButton.Enable = 'off';
                app.EvaluateStressAlignmentScaleButton.Enable = 'off';            
        end

        % Callback function
        function VoxelizeMeshGraphEdgesButtonPushed(app, event)
            global volumeFraction_;
            app.DataAlignmentButton.Enable = 'off';
            app.VoxelizeMeshGraphEdgesButton.Enable = 'off';
            pause(1);
            
            volumeFraction_ = MGD_VoxelizeMeshEdges();

            app.DesignVolumeFractionEditField.Value = volumeFraction_;

            app.DataAlignmentButton.Enable = 'on';
            app.VoxelizeMeshGraphEdgesButton.Enable = 'on';            
            app.GenerateVoxelbasedStructuralDesignButton.Enable = 'on';
            app.MainApp.ShowDesignbyDensityFieldMenu.Enable = 'on';
        end

        % Button pushed function: GenerateVoxelbasedStructuralDesignButton
        function GenerateVoxelbasedStructuralDesignButtonPushed(app, event)
            global volumeFractionDesign_;
            global outPath_;
            
            app.PreProcessPanel.Enable = 'off';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'off';
            app.GenerationSimulationPanel.Enable = 'off';
            app.DesignComplianceEditField.Value = 0;
            app.MainApp.DesignComplianceEditField.Value = 0;    
            pause(1);

            edgeWidth = app.EdgeThicknessEditField.Value; edgeWidth = max(edgeWidth,1);
            numLayerboundary = app.EntireBoundaryEditField.Value; numLayerboundary = max(numLayerboundary,0);
            numLayerLoads = app.LoadingAreaEditField.Value; numLayerLoads = max(numLayerLoads,0);
            numLayerFixation = app.FixedAreaEditField.Value; numLayerFixation = max(numLayerFixation,0);
            MGD_ConvertMeshGraph2MaterialLayout(edgeWidth, numLayerboundary, numLayerLoads, numLayerFixation);

            app.PreProcessPanel.Enable = 'on';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';
                app.GenerateVoxelbasedStructuralDesignButton.Enable = 'on';
                app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Enable = 'on';
                app.RestartLinearSystemSolvingButton.Enable = 'off';
                app.StressAnalysisonDesignButton.Enable = 'off';
                app.EvaluateStressAlignmentScaleButton.Enable = 'off';
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
            
            app.PreProcessPanel.Enable = 'off';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'off';
            app.GenerationSimulationPanel.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'off';    
            pause(1);

            GatherLSSandMPsettings(app.MainApp);
            [complianceDesign_, volumeFraction_] = FEA_ComputeComplianceVoxel(densityLayout_);
            
            app.PreProcessPanel.Enable = 'on';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';
                app.RestartLinearSystemSolvingButton.Enable = 'on';
                app.StressAnalysisonDesignButton.Enable = 'on';        
            app.ResultDisplayPanel.Enable = 'on'; 
                app.DesignComplianceEditField.Value = complianceDesign_;
            app.MainApp.DesignComplianceEditField.Value = complianceDesign_;
            ShowDeformation_Public(app.MainApp);          
        end

        % Button pushed function: EvaluateStressAlignmentScaleButton
        function EvaluateStressAlignmentScaleButtonPushed(app, event)
            global outPath_;
            global cartesianStressField_;
            
            disp('Compute Stress Aligment Scale between Solid and Design...');
            tStressAligmentAna = tic;
            dominantDirDesign = niftiread(strcat(outPath_, 'dominantDirDesign.nii'));                                                           
            alignmentMetricVolumeByEdgeAlignment = Common_ComputeEdgeAlignmentDeviation(dominantDirDesign);
            niftiwrite(alignmentMetricVolumeByEdgeAlignment, strcat(outPath_, 'alignmentMetricVolume_byEdge.nii'));            
            if ~isempty(cartesianStressField_)
                dominantDirSolid = niftiread(strcat(outPath_, 'dominantDirSolid.nii'));
                alignmentMetricVolumeByStressAlignment = Common_ComputeStressAlignmentDeviation(dominantDirSolid, dominantDirDesign);
                niftiwrite(alignmentMetricVolumeByStressAlignment, strcat(outPath_, 'alignmentMetricVolume_byStress.nii')); 
            end            
            disp(['Done with Stress Alignment Analysis after ', sprintf('%.1f', toc(tStressAligmentAna)), 's']);
            %system('"./src/quokka.exe" ./out/alignmentMetricVolume_byStress.nii'); 
        end

        % Button pushed function: RestartLinearSystemSolvingButton
        function RestartLinearSystemSolvingButtonPushed(app, event)
            global meshHierarchy_;
            global complianceDesign_;

            app.PreProcessPanel.Enable = 'off';
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
            
            app.PreProcessPanel.Enable = 'on';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';              
            app.ResultDisplayPanel.Enable = 'on'; 
                app.DesignComplianceEditField.Value = complianceDesign_;
            app.MainApp.DesignComplianceEditField.Value = complianceDesign_;        
        end

        % Button pushed function: StressAnalysisonDesignButton
        function StressAnalysisonDesignButtonPushed(app, event)
            global outPath_;            

            app.PreProcessPanel.Enable = 'off';
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
            
            app.PreProcessPanel.Enable = 'on';
            app.SettingsforMaterialLayoutConvertionPanel.Enable = 'on';
            app.GenerationSimulationPanel.Enable = 'on';
                app.EvaluateStressAlignmentScaleButton.Enable = 'on';
            app.ResultDisplayPanel.Enable = 'on';
        end

        % Menu selected function: SolidHexTetMeshmeshvtkmshMenu
        function SolidHexTetMeshmeshvtkmshMenuSelected(app, event)
            global axHandle_;
            global inputSolidMesh_;
            global frameStruct4Voxelization_;

            [fileName, dataPath] = uigetfile({'*.mesh'; '*.msh'; '*.vtk'}, 'Select a Solid Mesh File to Open');
            if isnumeric(fileName) || isnumeric(dataPath), return; end
            [~,~,fileExtension] = fileparts(fileName);
            if ~(strcmp(fileExtension, '.vtk') || strcmp(fileExtension, '.mesh') || strcmp(fileExtension, '.msh') || strcmp(fileExtension, '.obj'))
                warning('Un-supported Mesh/Graph Format!');
                return;
            end
            inputfileName = strcat(dataPath,fileName);
            IO_ImportSolidMesh(inputfileName);
            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            [az, el] = view(axHandle_);
            cla(axHandle_); colorbar(axHandle_, 'off');
            Vis_DrawMesh3D(axHandle_, inputSolidMesh_.boundaryNodeCoords, inputSolidMesh_.boundaryPatchNodMat, 1);
            view(axHandle_, az, el);

            frameStruct4Voxelization_ = Common_ExtractEdgeGraphFromSolidMeshTetHex(inputSolidMesh_);
            app.PreProcessPanel.Enable = 'on';
        end

        % Menu selected function: GraphobjMenu
        function GraphobjMenuSelected(app, event)
            global axHandle_;
            global vertexEdgeGraph_;
            global frameStruct4Voxelization_;

            [fileName, dataPath] = uigetfile({'*.obj'}, 'Select a Graph File to Open');
            if isnumeric(fileName) || isnumeric(dataPath), return; end
            [~,~,fileExtension] = fileparts(fileName);
            if ~strcmp(fileExtension, '.obj')
                warning('Un-supported Mesh/Graph Format!');
                return;
            end
            inputfileName = strcat(dataPath,fileName);
            IO_ImportVertexEdgeGraph(inputfileName);                
            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            [az, el] = view(axHandle_);
            cla(axHandle_); colorbar(axHandle_, 'off');
            Vis_DrawGraph3D(axHandle_, vertexEdgeGraph_.nodeCoords, vertexEdgeGraph_.eNodMat);
            view(axHandle_, az, el);
            Vis_UserLighting(axHandle_);
            frameStruct4Voxelization_ = vertexEdgeGraph_;

            app.PreProcessPanel.Enable = 'on';
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 649 591];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create OpenMenu
            app.OpenMenu = uimenu(app.UIFigure);
            app.OpenMenu.Text = 'Open';

            % Create SolidHexTetMeshmeshvtkmshMenu
            app.SolidHexTetMeshmeshvtkmshMenu = uimenu(app.OpenMenu);
            app.SolidHexTetMeshmeshvtkmshMenu.MenuSelectedFcn = createCallbackFcn(app, @SolidHexTetMeshmeshvtkmshMenuSelected, true);
            app.SolidHexTetMeshmeshvtkmshMenu.Text = 'Solid Hex/Tet Mesh ("*.mesh", "*.vtk", "*.msh")';

            % Create GraphobjMenu
            app.GraphobjMenu = uimenu(app.OpenMenu);
            app.GraphobjMenu.MenuSelectedFcn = createCallbackFcn(app, @GraphobjMenuSelected, true);
            app.GraphobjMenu.Text = 'Graph ("*.obj")';

            % Create EvaluateExternalMeshGraphbasedStructuralDesignPanel
            app.EvaluateExternalMeshGraphbasedStructuralDesignPanel = uipanel(app.UIFigure);
            app.EvaluateExternalMeshGraphbasedStructuralDesignPanel.Title = 'Evaluate External Mesh/Graph-based Structural Design';
            app.EvaluateExternalMeshGraphbasedStructuralDesignPanel.Position = [0 12 640 580];

            % Create SettingsforMaterialLayoutConvertionPanel
            app.SettingsforMaterialLayoutConvertionPanel = uipanel(app.EvaluateExternalMeshGraphbasedStructuralDesignPanel);
            app.SettingsforMaterialLayoutConvertionPanel.Title = 'Settings for Material Layout Convertion';
            app.SettingsforMaterialLayoutConvertionPanel.Position = [0 362 640 110];

            % Create EdgeThicknessEditFieldLabel
            app.EdgeThicknessEditFieldLabel = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.EdgeThicknessEditFieldLabel.HorizontalAlignment = 'right';
            app.EdgeThicknessEditFieldLabel.Position = [445 16 90 22];
            app.EdgeThicknessEditFieldLabel.Text = 'Edge Thickness';

            % Create EdgeThicknessEditField
            app.EdgeThicknessEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.EdgeThicknessEditField.ValueDisplayFormat = '%.0f';
            app.EdgeThicknessEditField.Position = [550 16 68 22];
            app.EdgeThicknessEditField.Value = 3;

            % Create EntireBoundaryEditFieldLabel
            app.EntireBoundaryEditFieldLabel = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.EntireBoundaryEditFieldLabel.HorizontalAlignment = 'right';
            app.EntireBoundaryEditFieldLabel.Position = [444 56 91 22];
            app.EntireBoundaryEditFieldLabel.Text = 'Entire Boundary';

            % Create EntireBoundaryEditField
            app.EntireBoundaryEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.EntireBoundaryEditField.ValueDisplayFormat = '%.0f';
            app.EntireBoundaryEditField.Position = [550 56 68 22];
            app.EntireBoundaryEditField.Value = 2;

            % Create LoadingAreaEditField_4Label
            app.LoadingAreaEditField_4Label = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.LoadingAreaEditField_4Label.HorizontalAlignment = 'right';
            app.LoadingAreaEditField_4Label.Position = [16 54 76 22];
            app.LoadingAreaEditField_4Label.Text = 'Loading Area';

            % Create LoadingAreaEditField
            app.LoadingAreaEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.LoadingAreaEditField.ValueDisplayFormat = '%.0f';
            app.LoadingAreaEditField.Position = [107 54 68 22];

            % Create FixedAreaEditField_4Label
            app.FixedAreaEditField_4Label = uilabel(app.SettingsforMaterialLayoutConvertionPanel);
            app.FixedAreaEditField_4Label.HorizontalAlignment = 'right';
            app.FixedAreaEditField_4Label.Position = [240 55 62 22];
            app.FixedAreaEditField_4Label.Text = 'Fixed Area';

            % Create FixedAreaEditField
            app.FixedAreaEditField = uieditfield(app.SettingsforMaterialLayoutConvertionPanel, 'numeric');
            app.FixedAreaEditField.ValueDisplayFormat = '%.0f';
            app.FixedAreaEditField.Position = [317 55 68 22];

            % Create GenerationSimulationPanel
            app.GenerationSimulationPanel = uipanel(app.EvaluateExternalMeshGraphbasedStructuralDesignPanel);
            app.GenerationSimulationPanel.Title = 'Generation & Simulation';
            app.GenerationSimulationPanel.Position = [0 145 640 218];

            % Create StiffnessEvaluationofVoxelbasedStructuralDesignButton
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.ButtonPushedFcn = createCallbackFcn(app, @StiffnessEvaluationofVoxelbasedStructuralDesignButtonPushed, true);
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Position = [320 113 302 23];
            app.StiffnessEvaluationofVoxelbasedStructuralDesignButton.Text = 'Stiffness Evaluation of Voxel based Structural Design';

            % Create EvaluateStressAlignmentScaleButton
            app.EvaluateStressAlignmentScaleButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.EvaluateStressAlignmentScaleButton.ButtonPushedFcn = createCallbackFcn(app, @EvaluateStressAlignmentScaleButtonPushed, true);
            app.EvaluateStressAlignmentScaleButton.Position = [434 17 188 23];
            app.EvaluateStressAlignmentScaleButton.Text = 'Evaluate Stress Alignment Scale';

            % Create RestartLinearSystemSolvingButton
            app.RestartLinearSystemSolvingButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.RestartLinearSystemSolvingButton.ButtonPushedFcn = createCallbackFcn(app, @RestartLinearSystemSolvingButtonPushed, true);
            app.RestartLinearSystemSolvingButton.Position = [54 113 180 23];
            app.RestartLinearSystemSolvingButton.Text = 'Re-start Linear System Solving';

            % Create StressAnalysisonDesignButton
            app.StressAnalysisonDesignButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.StressAnalysisonDesignButton.ButtonPushedFcn = createCallbackFcn(app, @StressAnalysisonDesignButtonPushed, true);
            app.StressAnalysisonDesignButton.Position = [468 66 154 23];
            app.StressAnalysisonDesignButton.Text = 'Stress Analysis on Design';

            % Create GenerateVoxelbasedStructuralDesignButton
            app.GenerateVoxelbasedStructuralDesignButton = uibutton(app.GenerationSimulationPanel, 'push');
            app.GenerateVoxelbasedStructuralDesignButton.ButtonPushedFcn = createCallbackFcn(app, @GenerateVoxelbasedStructuralDesignButtonPushed, true);
            app.GenerateVoxelbasedStructuralDesignButton.BackgroundColor = [0.9608 0.9608 0.9608];
            app.GenerateVoxelbasedStructuralDesignButton.Position = [392 156 231 23];
            app.GenerateVoxelbasedStructuralDesignButton.Text = 'Generate Voxel based Structural Design';

            % Create PreProcessPanel
            app.PreProcessPanel = uipanel(app.EvaluateExternalMeshGraphbasedStructuralDesignPanel);
            app.PreProcessPanel.Title = 'Pre-Process';
            app.PreProcessPanel.Position = [0 471 640 85];

            % Create DataAlignmentButton
            app.DataAlignmentButton = uibutton(app.PreProcessPanel, 'push');
            app.DataAlignmentButton.ButtonPushedFcn = createCallbackFcn(app, @DataAlignmentButtonPushed, true);
            app.DataAlignmentButton.Position = [519 22 100 23];
            app.DataAlignmentButton.Text = 'Data Alignment';

            % Create ResultDisplayPanel
            app.ResultDisplayPanel = uipanel(app.EvaluateExternalMeshGraphbasedStructuralDesignPanel);
            app.ResultDisplayPanel.Title = 'Result Display';
            app.ResultDisplayPanel.Position = [0 56 640 90];

            % Create DesignComplianceEditFieldLabel
            app.DesignComplianceEditFieldLabel = uilabel(app.ResultDisplayPanel);
            app.DesignComplianceEditFieldLabel.HorizontalAlignment = 'right';
            app.DesignComplianceEditFieldLabel.Position = [28 19 109 22];
            app.DesignComplianceEditFieldLabel.Text = 'Design Compliance';

            % Create DesignComplianceEditField
            app.DesignComplianceEditField = uieditfield(app.ResultDisplayPanel, 'numeric');
            app.DesignComplianceEditField.Position = [152 19 100 22];

            % Create DesignVolumeFractionEditFieldLabel
            app.DesignVolumeFractionEditFieldLabel = uilabel(app.ResultDisplayPanel);
            app.DesignVolumeFractionEditFieldLabel.HorizontalAlignment = 'right';
            app.DesignVolumeFractionEditFieldLabel.Position = [370 19 132 22];
            app.DesignVolumeFractionEditFieldLabel.Text = 'Design Volume Fraction';

            % Create DesignVolumeFractionEditField
            app.DesignVolumeFractionEditField = uieditfield(app.ResultDisplayPanel, 'numeric');
            app.DesignVolumeFractionEditField.Position = [517 19 100 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Mtd_MeshGraphBasedStructDesign(varargin)

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