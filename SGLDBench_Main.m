classdef SGLDBench_Main < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        FileMenu                        matlab.ui.container.Menu
        ImportMenu                      matlab.ui.container.Menu
        TriangularSurfaceMeshplyobjstlMenu  matlab.ui.container.Menu
        VoxelModelTopVoxelMenu          matlab.ui.container.Menu
        BuiltinShapesMenu               matlab.ui.container.Menu
        CuboidMenu                      matlab.ui.container.Menu
        LshapeMenu                      matlab.ui.container.Menu
        CylinderMenu                    matlab.ui.container.Menu
        ExportMenu                      matlab.ui.container.Menu
        DesignVolumeniiMenu             matlab.ui.container.Menu
        VoxelModelTopVoxelMenu_2        matlab.ui.container.Menu
        VoxelModelwithDensityLayoutTopVoxelMenu  matlab.ui.container.Menu
        VisualizationMenu               matlab.ui.container.Menu
        ShowDesignVolumeWebGLMenu       matlab.ui.container.Menu
        ShowInputTriangularSurfaceMeshMenu  matlab.ui.container.Menu
        ShowProblemDescriptionMenu      matlab.ui.container.Menu
        ShowDesignDomainMenu            matlab.ui.container.Menu
        ShowDeformationMenu             matlab.ui.container.Menu
        ShowStressFieldvonMisesStressMenu  matlab.ui.container.Menu
        ShowPSLsMenu                    matlab.ui.container.Menu
        ShowVertexEdgeGraphMenu         matlab.ui.container.Menu
        ShowDesignbyIsosurfaceNotrecommendedMenu  matlab.ui.container.Menu
        TabGroup3                       matlab.ui.container.TabGroup
        ModelingTab                     matlab.ui.container.Tab
        ApplyforBoundaryConditionsPanel  matlab.ui.container.Panel
        TabGroupSelection               matlab.ui.container.TabGroup
        BoxSelectionTab                 matlab.ui.container.Tab
        CornerPot2ZEditField            matlab.ui.control.NumericEditField
        CornerPot2ZEditFieldLabel       matlab.ui.control.Label
        CornerPot1ZEditField            matlab.ui.control.NumericEditField
        CornerPot1ZEditFieldLabel       matlab.ui.control.Label
        CornerPot2YEditField            matlab.ui.control.NumericEditField
        CornerPot2YLabel                matlab.ui.control.Label
        CornerPot1YEditField            matlab.ui.control.NumericEditField
        CornerPot1YEditFieldLabel       matlab.ui.control.Label
        CornerPot2XEditField            matlab.ui.control.NumericEditField
        CornerPot2XLabel                matlab.ui.control.Label
        CornerPot1XEditField            matlab.ui.control.NumericEditField
        CornerPot1XLabel                matlab.ui.control.Label
        SphereSelectionTab              matlab.ui.container.Tab
        RadiusEditField                 matlab.ui.control.NumericEditField
        RadiusEditFieldLabel            matlab.ui.control.Label
        CenterZEditField                matlab.ui.control.NumericEditField
        CenterZEditFieldLabel           matlab.ui.control.Label
        CenterYEditField                matlab.ui.control.NumericEditField
        CenterYEditFieldLabel           matlab.ui.control.Label
        CenterXEditField                matlab.ui.control.NumericEditField
        CenterXEditFieldLabel           matlab.ui.control.Label
        SelectionOptionsDropDown        matlab.ui.control.DropDown
        SelectionOptionsDropDownLabel   matlab.ui.control.Label
        BuiltinBoundaryConditionsPanel  matlab.ui.container.Panel
        OnlyCuboidDomainDropDown        matlab.ui.control.DropDown
        OnlyCuboidDomainDropDownLabel   matlab.ui.control.Label
        TabGroup2                       matlab.ui.container.TabGroup
        FixingTab                       matlab.ui.container.Tab
        XDirFixedCheckBox               matlab.ui.control.CheckBox
        YDirFixedCheckBox               matlab.ui.control.CheckBox
        ClearFixationButton             matlab.ui.control.Button
        ApplyforFixationButton          matlab.ui.control.Button
        ZDirFixedCheckBox               matlab.ui.control.CheckBox
        LoadingTab                      matlab.ui.container.Tab
        FzNEditField                    matlab.ui.control.NumericEditField
        FzNEditFieldLabel               matlab.ui.control.Label
        FyNEditField                    matlab.ui.control.NumericEditField
        FyNEditFieldLabel               matlab.ui.control.Label
        ClearLoadsButton                matlab.ui.control.Button
        ApplyforLoadsButton             matlab.ui.control.Button
        FxNEditField                    matlab.ui.control.NumericEditField
        FxLabel                         matlab.ui.control.Label
        NodeSelectionButton_2           matlab.ui.control.Button
        NodeUnSelectionButton           matlab.ui.control.Button
        DomainVoxelizationPanel         matlab.ui.container.Panel
        DOFsEditField                   matlab.ui.control.NumericEditField
        DOFsEditFieldLabel              matlab.ui.control.Label
        ElementsEditField               matlab.ui.control.NumericEditField
        ElementsEditFieldLabel          matlab.ui.control.Label
        VoxelizingButton                matlab.ui.control.Button
        CellSizeEditField               matlab.ui.control.NumericEditField
        CellSizeEditFieldLabel          matlab.ui.control.Label
        TargetVoxelResolutionEditField  matlab.ui.control.NumericEditField
        TargetVoxelResolutionEditFieldLabel  matlab.ui.control.Label
        SimulationTab                   matlab.ui.container.Tab
        AdditionalStiffnessEvaluationPanel  matlab.ui.container.Panel
        StiffnesswrtChangedLoadingDirectionsButton  matlab.ui.control.Button
        StiffnessofExternalMeshGraphviaVoxelFEAButton  matlab.ui.control.Button
        StressSimulationonSolidDomainPanel  matlab.ui.container.Panel
        DesignVolEditField              matlab.ui.control.NumericEditField
        DesignVolEditFieldLabel         matlab.ui.control.Label
        DesignComplianceEditField       matlab.ui.control.NumericEditField
        DesignComplianceEditFieldLabel  matlab.ui.control.Label
        FEAwithSolidDesignDomainButton  matlab.ui.control.Button
        SolidComplianceEditField        matlab.ui.control.NumericEditField
        SolidComplianceEditFieldLabel   matlab.ui.control.Label
        MaterialLayoutDesignPanel       matlab.ui.container.Panel
        SimulationTasksDropDown         matlab.ui.control.DropDown
        SimulationTasksDropDownLabel    matlab.ui.control.Label
        LinearSystemSolverPanel         matlab.ui.container.Panel
        WeightingFactorofJacobiSmoothingProcessEditField  matlab.ui.control.NumericEditField
        WeightingFactorofJacobiSmoothingProcessEditFieldLabel  matlab.ui.control.Label
        MaximumIterationsEditField      matlab.ui.control.NumericEditField
        MaximumIterationsEditFieldLabel  matlab.ui.control.Label
        ResidualThresholdEditField      matlab.ui.control.NumericEditField
        ResidualThresholdEditFieldLabel  matlab.ui.control.Label
        MaterialPropertiesPanel         matlab.ui.container.Panel
        YoungsModulusCompliantmaterialEditField  matlab.ui.control.NumericEditField
        YoungsModulusCompliantmaterialEditFieldLabel  matlab.ui.control.Label
        PoissonsRatioEditField          matlab.ui.control.NumericEditField
        PoissonsRatioEditFieldLabel     matlab.ui.control.Label
        YoungsModulusStiffmaterialEditField  matlab.ui.control.NumericEditField
        YoungsModulusStiffmaterialEditFieldLabel  matlab.ui.control.Label
    end

    
    properties (Access = private)
        %%app_LinearSystemSolver_settings % Reference to LSS settings
        app_ObjectiSelectionWindow_settings % Reference to Object Selection
        %%app_Optimizer_settings % Reference to Optimizer settings
        %%app_MaterialProperties_settings % Reference to Material Property settings
        comp_SimTask_TopologyOptimization_Func
        comp_SimTask_PSLsGuidedStructDesign_Func
        comp_SimTask_MeshGraphBasedStructDesign_Func
        comp_SimTask_DesignWRTchangedLoadingDirections_Func
        comp_SimTask_StressAwareGeoSyn3rdParty_AroraMethod
        comp_SimTask_StressAwareGeoSyn3rdParty_GaoMethod
        comp_SimTask_StressAwareGeoSyn3rdParty_VoronoiMethod
        app_Mdl_CreateCuboidDesignDomain
        app_Mdl_CreateLshapeDesignDomain
        app_Mdl_CreateSphereShellDesignDomain
        app_Mdl_CreateCylinderDesignDomain
    end
    
    methods (Access = public)
        
        function ShowVoxelizedModel_Public(app)
            InputVoxelsMenuSelected(app);
        end

        function ShowProblemDescription_Public(app)
            ShowProblemDescriptionMenuSelected(app)
        end

        function ShowTriSurfaceMesh_public(app)
            ShowInputTriSurfaceMeshMenuSelected(app);
        end

        function ShowDeformation_Public(app)
            TotalMenu_2Selected(app);
        end

        function ShowStressComponents_Public(app)
            VonMisesMenuSelected(app);
        end

        function ShowPSLs_Public(app)
            PSLsMenuSelected(app);
        end

        function ShowDesignDomain_Public(app)
            ShowDesignDomainMenuSelected(app);
        end
        
        function ShowDesignDensityLayout_Public(app)
            ShowDesignbyIsosurfaceNotrecommendedMenuSelected(app)
        end

        function FEA_Public(app)
            FEAwithSolidDesignDomainButtonPushed(app)
        end

        function StressAnalysis_Public(app)
            StressAnalysisButtonPushed(app)
        end
        
        function ShowEdgeVertexGraph_Public(app)
            InputEdgeVertexGraphMenuSelected(app)
        end

        function SetupSelectionOptions_Public(app)
            global boundingBox_;
            app.CornerPot1XEditField.Value = boundingBox_(1,1); 
            app.CornerPot1YEditField.Value = boundingBox_(1,2); 
            app.CornerPot1ZEditField.Value = boundingBox_(1,3);
            app.CornerPot2XEditField.Value = boundingBox_(2,1); 
            app.CornerPot2YEditField.Value = boundingBox_(2,2); 
            app.CornerPot2ZEditField.Value = boundingBox_(2,3);
            
            app.CenterXEditField.Value = boundingBox_(2,1);
            app.CenterYEditField.Value = boundingBox_(2,2);
            app.CenterZEditField.Value = boundingBox_(2,3);
            app.RadiusEditField.Value = 10;
        end        
        
        function MainWindowCtrl(app, opt)
            if opt
                app.FileMenu.Enable = 'on';
                app.VisualizationMenu.Enable = 'on';
                app.DomainVoxelizationPanel.Enable = 'on';
                app.ApplyforBoundaryConditionsPanel.Enable = 'on';
                app.MaterialPropertiesPanel.Enable = 'on';
                app.StressSimulationonSolidDomainPanel.Enable = 'on';
                app.MaterialLayoutDesignPanel.Enable = 'on';
                app.AdditionalStiffnessEvaluationPanel.Enable = 'on';                    
            else
                app.FileMenu.Enable = 'off';
                app.VisualizationMenu.Enable = 'off';
                app.DomainVoxelizationPanel.Enable = 'off';
                app.ApplyforBoundaryConditionsPanel.Enable = 'off';
                app.MaterialPropertiesPanel.Enable = 'off';
                app.StressSimulationonSolidDomainPanel.Enable = 'off';
                app.MaterialLayoutDesignPanel.Enable = 'off';
                app.AdditionalStiffnessEvaluationPanel.Enable = 'off';
            end
        end

        function GatherLSSandMPsettings(app)
            global tol_;
            global maxIT_;
            global weightFactorJacobi_;
            % global nonDyadic_;
            % global MEXfunc_;
            global modulus_;
            global poissonRatio_;
            global modulusMin_;
            global cellSize_;

            tol_ = app.ResidualThresholdEditField.Value;
            maxIT_ = app.MaximumIterationsEditField.Value;
            weightFactorJacobi_ = app.WeightingFactorofJacobiSmoothingProcessEditField.Value;
            % nonDyadic_ = app.NonDyadicCheckBox.Value;
            % MEXfunc_ = app.MEXFuncCheckBox.Value;
            modulus_ = app.YoungsModulusStiffmaterialEditField.Value;
            poissonRatio_ = app.PoissonsRatioEditField.Value;            
            modulusMin_ = app.YoungsModulusCompliantmaterialEditField.Value;
            cellSize_ = app.CellSizeEditField.Value;            
        end

        function InitializeAppParameters(app)
            global finestResolutionControl_;
            global cellSize_;
            global modulus_;
            global poissonRatio_;
            global modulusMin_;
            global tol_;
            global maxIT_;
            global weightFactorJacobi_;
            
            app.TargetVoxelResolutionEditField.Value = finestResolutionControl_;
            app.CellSizeEditField.Value = cellSize_;
            app.ElementsEditField.Value = 0;
            app.DOFsEditField.Value = 0;
            app.FxNEditField.Value = 0;
            app.FyNEditField.Value = 0;
            app.FzNEditField.Value = 0;
            app.OnlyCuboidDomainDropDown.Value = 'None';
            app.XDirFixedCheckBox.Value = 1;
            app.YDirFixedCheckBox.Value = 1;
            app.ZDirFixedCheckBox.Value = 1;
            app.YoungsModulusCompliantmaterialEditField.Value = modulusMin_;
            app.YoungsModulusStiffmaterialEditField.Value = modulus_;
            app.PoissonsRatioEditField.Value = poissonRatio_;
            app.ResidualThresholdEditField.Value = tol_;
            app.MaximumIterationsEditField.Value = maxIT_;
            app.WeightingFactorofJacobiSmoothingProcessEditField.Value = weightFactorJacobi_;
            app.SolidComplianceEditField.Value = 0;
            app.SelectionOptionsDropDown.Value = 'None';
        end

        function InitializeMainAppInterface(app)
            %%Prevent Accidental Touches            
            app.ExportMenu.Enable = 'off';
            app.DomainVoxelizationPanel.Enable = 'off';
                app.TargetVoxelResolutionEditField.Enable = 'on';
                app.CellSizeEditField.Enable = 'on';
                app.VoxelizingButton.Enable = 'on';
            app.ApplyforBoundaryConditionsPanel.Enable = 'off';                
            app.BuiltinBoundaryConditionsPanel.Enable = 'off';
            app.StressSimulationonSolidDomainPanel.Enable = 'off';
            app.MaterialLayoutDesignPanel.Enable = 'off';
            app.AdditionalStiffnessEvaluationPanel.Enable = 'off';
                app.StiffnessofExternalMeshGraphviaVoxelFEAButton.Enable = 'off';
                app.StiffnesswrtChangedLoadingDirectionsButton.Enable = 'off';
            app.VisualizationMenu.Enable = 'on';
                app.ShowInputTriangularSurfaceMeshMenu.Enable = 'off';
                app.ShowProblemDescriptionMenu.Enable = 'off';
                app.ShowDesignDomainMenu.Enable = 'off';
                app.ShowDeformationMenu.Enable = 'off';
                app.ShowStressFieldvonMisesStressMenu.Enable = 'off';
                app.ShowDesignbyIsosurfaceNotrecommendedMenu.Enable = 'off';
                app.ShowVertexEdgeGraphMenu.Enable = 'off';
                app.ShowPSLsMenu.Enable = 'off';
            % DisableSelectionTab(app);
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            clc;
            global axHandle_;            
            addpath('./GUIs/');
            addpath('./src/');
            addpath('./src/MEXfuncs/');
            if ~exist('./out/', 'dir'), mkdir('./out/'); end
            %%Initialization
            Data_GlobalVariables;
            outPath_ = './out/';
            InitializeMainAppInterface(app);
            InitializeAppParameters(app);
            close all;
            figure; axHandle_ = gca; view(axHandle_,3);
            cla(axHandle_);
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            delete(app.app_ObjectiSelectionWindow_settings)
            delete(app.comp_SimTask_DesignWRTchangedLoadingDirections_Func)
            delete(app.comp_SimTask_TopologyOptimization_Func)
            delete(app.comp_SimTask_PSLsGuidedStructDesign_Func)
            delete(app.comp_SimTask_MeshGraphBasedStructDesign_Func)
            delete(app.comp_SimTask_StressAwareGeoSyn3rdParty_AroraMethod)
            delete(app.comp_SimTask_StressAwareGeoSyn3rdParty_GaoMethod)
            delete(app.comp_SimTask_StressAwareGeoSyn3rdParty_VoronoiMethod)
            delete(app.app_Mdl_CreateCuboidDesignDomain)
            delete(app.app_Mdl_CreateLshapeDesignDomain)
            delete(app.app_Mdl_CreateCylinderDesignDomain)
            delete(app)            
        end

        % Menu selected function: CuboidMenu
        function CuboidMenuSelected(app, event)
            global axHandle_;
            %%Reset App            
            Data_GlobalVariables;
            InitializeMainAppInterface(app);
            InitializeAppParameters(app);  

            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            cla(axHandle_);
            app.app_Mdl_CreateCuboidDesignDomain = Mdl_Cuboid(app);
            
            app.ImportMenu.Enable = 'off';
            % app.TargetVoxelResolutionEditField.Enable = 'on';
            app.ExportMenu.Enable = 'off';
            % app.VisualizationMenu.Enable = 'on';
            % app.DomainVoxelizationPanel.Enable = 'on';
            app.ApplyforBoundaryConditionsPanel.Enable = 'off';
                app.BuiltinBoundaryConditionsPanel.Enable = 'on';
            app.StressSimulationonSolidDomainPanel.Enable = 'off';
            app.MaterialLayoutDesignPanel.Enable = 'off';
        end

        % Value changed function: OnlyCuboidDomainDropDown
        function OnlyCuboidDomainDropDownValueChanged(app, event)
            value = app.OnlyCuboidDomainDropDown.Value;
            FEA_BuiltInBoundaryConditions4CuboidDesignDomain(value);
            ShowProblemDescriptionMenuSelected(app, event);
            % app.SimulationTasksDropDown.Enable = 'on';
            % app.EvaluationTasksDropDown.Enable = 'on';
            % app.FEAwithSolidDesignDomainButton.Enable = 'on';
            app.StressSimulationonSolidDomainPanel.Enable = 'on';
            app.MaterialLayoutDesignPanel.Enable = 'on';
            app.AdditionalStiffnessEvaluationPanel.Enable = 'on';
                app.StiffnessofExternalMeshGraphviaVoxelFEAButton.Enable = 'on';
        end

        % Menu selected function: LshapeMenu
        function LshapeMenuSelected(app, event)
            global axHandle_;
            %%Reset App            
            Data_GlobalVariables;
            InitializeMainAppInterface(app);
            InitializeAppParameters(app);  

            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            cla(axHandle_);
            app.app_Mdl_CreateLshapeDesignDomain = Mdl_Lshape(app);
            
            app.ImportMenu.Enable = 'off';
            % app.TargetVoxelResolutionEditField.Enable = 'on';
            app.ExportMenu.Enable = 'off';
            % app.VisualizationMenu.Enable = 'on';
            % app.DomainVoxelizationPanel.Enable = 'on';
            app.ApplyforBoundaryConditionsPanel.Enable = 'off';
            app.StressSimulationonSolidDomainPanel.Enable = 'off';
            app.MaterialLayoutDesignPanel.Enable = 'off';           
        end

        % Menu selected function: CylinderMenu
        function CylinderMenuSelected(app, event)
            global axHandle_;
            %%Reset App            
            Data_GlobalVariables;
            InitializeMainAppInterface(app);
            InitializeAppParameters(app);  
           
            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            cla(axHandle_);
            app.app_Mdl_CreateCylinderDesignDomain = Mdl_Cylinder(app);
            
            app.ImportMenu.Enable = 'off';
            % app.TargetVoxelResolutionEditField.Enable = 'on';
            app.ExportMenu.Enable = 'off';
            % app.VisualizationMenu.Enable = 'on';
            % app.DomainVoxelizationPanel.Enable = 'on';
            app.ApplyforBoundaryConditionsPanel.Enable = 'off';
            app.StressSimulationonSolidDomainPanel.Enable = 'off';
            app.MaterialLayoutDesignPanel.Enable = 'off';          
        end

        % Menu selected function: TriangularSurfaceMeshplyobjstlMenu
        function SurfaceTriMeshplyobjMenuSelected(app, event)
            global axHandle_;
            %%Reset App
            Data_GlobalVariables;
            InitializeMainAppInterface(app);
            InitializeAppParameters(app);

            [fileName, dataPath] = uigetfile({'*.ply'; '*.obj'; '*.stl'}, 'Select a Surface Mesh File to Open');
            if isnumeric(fileName) || isnumeric(dataPath), return; end
            [~,~,fileExtension] = fileparts(fileName);
            if ~(strcmp(fileExtension, '.ply') || strcmp(fileExtension, '.obj') || strcmp(fileExtension, '.stl'))
                warning('Un-supported Mesh Format!');
                return;
            end
            inputSurfaceMeshfileName = strcat(dataPath,fileName);
            IO_ImportSurfaceMesh(inputSurfaceMeshfileName);

            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            cla(axHandle_);
            ShowInputTriSurfaceMeshMenuSelected(app, event);
            app.ShowInputTriangularSurfaceMeshMenu.Enable = 'on';
            app.DomainVoxelizationPanel.Enable = 'on';        

            app.TargetVoxelResolutionEditField.Editable = 'on';
            app.VoxelModelTopVoxelMenu.Enable = 'on';
            
        end

        % Menu selected function: VoxelModelTopVoxelMenu
        function ImportVoxelModelTopVoxelMenuSelected(app, event)
            global axHandle_;
            %%Reset App
            Data_GlobalVariables;
            InitializeMainAppInterface(app);
            InitializeAppParameters(app);  

            [fileName, dataPath] = uigetfile('*.TopVoxel', 'Select a Voxel File to Open');
            if isnumeric(fileName) || isnumeric(dataPath), return; end
            [~,~,fileExtension] = fileparts(fileName);
            if ~strcmp(fileExtension, '.TopVoxel')
                warning('Un-supported Mesh Format!');
                return;
            end
            inputVoxelfileName = strcat(dataPath,fileName);
            IO_ImportTopVoxels(inputVoxelfileName);
            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            cla(axHandle_);
            ShowProblemDescriptionMenuSelected(app, event);
            
            app.VisualizationMenu.Enable = 'on';
                app.ShowProblemDescriptionMenu.Enable = 'on';
                app.ShowDesignDomainMenu.Enable = 'on';
                app.ExportMenu.Enable = 'on';            
            app.DomainVoxelizationPanel.Enable = 'on';
                app.ElementsEditField.Value = meshHierarchy_(1).numElements;
                app.ElementsEditField.Editable = 'off';
                app.DOFsEditField.Value = meshHierarchy_(1).numDOFs;
                app.DOFsEditField.Editable = 'off';
                app.TargetVoxelResolutionEditField.Value = max([meshHierarchy_(1).resX meshHierarchy_(1).resY meshHierarchy_(1).resZ]);
                app.TargetVoxelResolutionEditField.Editable = 'off';
                app.VoxelizingButton.Enable = 'off';
            % EnableSelectionTab(app);
            app.ApplyforBoundaryConditionsPanel.Enable = 'on';
            
            %%Method 3, 5, 6 does not work for voxel input
            app.SimulationTasksDropDown.Items = {'None', 'Mtd - Topology Optimization', 'Mtd - Porous Infill Optimization', 'Mtd - PSLs-guided Infill Design'};       

            SetupSelectionOptions_Public(app);
            if ~(isempty(loadingCond_) || isempty(fixingCond_))
                app.StressSimulationonSolidDomainPanel.Enable = 'on';
                app.MaterialLayoutDesignPanel.Enable = 'on';
                app.AdditionalStiffnessEvaluationPanel.Enable = 'on';
                    app.StiffnessofExternalMeshGraphviaVoxelFEAButton.Enable = 'on';
            end
            if ~isempty(densityLayout_)
                app.ShowDesignbyIsosurfaceNotrecommendedMenu.Enable = 'on';
            end
        end

        % Button pushed function: VoxelizingButton
        function VoxelizingButtonPushed(app, event)
            global meshHierarchy_;
            % global nelx_; global nely_; global nelz_;
            % global boundingBox_;
            
            app.DomainVoxelizationPanel.Enable = 'off';
            app.FileMenu.Enable = 'off';
            app.VisualizationMenu.Enable = 'off';
            app.ApplyforBoundaryConditionsPanel.Enable = 'off';
            app.StressSimulationonSolidDomainPanel.Enable = 'off';
            app.MaterialLayoutDesignPanel.Enable = 'off';
            pause(1);

            tStart = tic;
            FEA_CreateVoxelizedModel(app.TargetVoxelResolutionEditField.Value);
            disp(['Voxelizing Domain Costs: ' sprintf('%10.3g',toc(tStart)) 's']);
            tStart = tic;
            FEA_VoxelBasedDiscretization();
            disp(['Setup Voxel-based FEA Model: ' sprintf('%10.3g',toc(tStart)) 's']);
            ShowProblemDescriptionMenuSelected(app, event);

            app.DomainVoxelizationPanel.Enable = 'on';
                app.ElementsEditField.Value = meshHierarchy_(1).numElements;
                app.DOFsEditField.Value = meshHierarchy_(1).numDOFs;          
            app.FileMenu.Enable = 'on';
                app.ExportMenu.Enable = 'on';
            app.VisualizationMenu.Enable = 'on';
                app.ShowProblemDescriptionMenu.Enable = 'on';
                app.ShowDesignDomainMenu.Enable = 'on';
            % EnableSelectionTab(app);    
            app.ApplyforBoundaryConditionsPanel.Enable = 'on';                
                SetupSelectionOptions_Public(app);
            app.SimulationTasksDropDown.Items = {'None', 'Mtd - Topology Optimization', 'Mtd - Porous Infill Optimization', ...
                'Mtd - Stress-aware Graded Voronoi Diagram Infill Design', 'Mtd - PSLs-guided Infill Design', ...
                    'Mtd - Stress-aligned Conforming Lattice Infill Design', 'Mtd - Stress-aligned Volumetric Michell Trusses Infill Design'};            
        end

        % Menu selected function: ShowInputTriangularSurfaceMeshMenu
        function ShowInputTriSurfaceMeshMenuSelected(app, event)
            global axHandle_;
            global surfaceTriMesh_;
            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            [az, el] = view(axHandle_);
            cla(axHandle_); colorbar(axHandle_, 'off'); 
            Vis_DrawMesh3D(axHandle_, surfaceTriMesh_.nodeCoords, surfaceTriMesh_.eNodMat, 1);
            view(axHandle_, az, el);
        end

        % Button pushed function: ApplyforLoadsButton
        function ApplyforLoadsButtonPushed(app, event)
            global axHandle_;
            global loadingCond_;
            global fixingCond_;             
            if 0==app.FxNEditField.Value && 0==app.FyNEditField.Value && 0==app.FzNEditField.Value, return; end
            forceVec = [app.FxNEditField.Value app.FyNEditField.Value app.FzNEditField.Value];
            iLoadingVec2Draw = FEA_Apply4Loads(forceVec);
            if isempty(iLoadingVec2Draw), return; end
            Interaction_ClearPickedNodes();
            Vis_ShowLoadingCondition(axHandle_, iLoadingVec2Draw);
            if ~isempty(loadingCond_) && ~isempty(fixingCond_)
                app.StressSimulationonSolidDomainPanel.Enable = 'on';
                app.MaterialLayoutDesignPanel.Enable = 'on';
                app.AdditionalStiffnessEvaluationPanel.Enable = 'on';
                    app.StiffnessofExternalMeshGraphviaVoxelFEAButton.Enable = 'on';
            end
        end

        % Button pushed function: ClearLoadsButton
        function ClearLoadsButtonPushed(app, event)
            global loadingCond_;
            global F_; F_ = [];
            loadingCond_ = [];
            app.StressSimulationonSolidDomainPanel.Enable = 'off';
            app.MaterialLayoutDesignPanel.Enable = 'off';
            app.AdditionalStiffnessEvaluationPanel.Enable = 'off';            
            ShowProblemDescriptionMenuSelected(app, event);
        end

        % Button pushed function: ApplyforFixationButton
        function ApplyforFixationButtonPushed(app, event)
            global axHandle_;
            global loadingCond_;
            global fixingCond_;           
            if ~(app.XDirFixedCheckBox.Value || app.YDirFixedCheckBox.Value || app.ZDirFixedCheckBox.Value), return; end
            fixingOpt = zeros(1,3);
            if app.XDirFixedCheckBox.Value, fixingOpt(1) = 1; end
            if app.YDirFixedCheckBox.Value, fixingOpt(2) = 1; end
            if app.ZDirFixedCheckBox.Value, fixingOpt(3) = 1; end
            iFixingArr2Draw = FEA_Apply4Fixations(fixingOpt);
            if isempty(iFixingArr2Draw), return; end
            Interaction_ClearPickedNodes();
            Vis_ShowFixingCondition(axHandle_, iFixingArr2Draw);
            if ~isempty(loadingCond_) && ~isempty(fixingCond_)
                app.StressSimulationonSolidDomainPanel.Enable = 'on';
                app.MaterialLayoutDesignPanel.Enable = 'on';
                app.AdditionalStiffnessEvaluationPanel.Enable = 'on';
                    app.StiffnessofExternalMeshGraphviaVoxelFEAButton.Enable = 'on';
            end            
        end

        % Menu selected function: ShowProblemDescriptionMenu
        function ShowProblemDescriptionMenuSelected(app, event)
            global axHandle_;
            global meshHierarchy_;
            global loadingCond_;
            global fixingCond_;
            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            [az, el] = view(axHandle_);
            cla(axHandle_); colorbar(axHandle_, 'off');
            Vis_DrawMesh3D(axHandle_, meshHierarchy_(1).boundaryNodeCoords, meshHierarchy_(1).boundaryEleFaces, 0);
            Vis_ShowLoadingCondition(axHandle_, loadingCond_);
            Vis_ShowFixingCondition(axHandle_, fixingCond_);
            view(axHandle_, az, el);
            axis(axHandle_, 'on'); xlabel('X'); ylabel('Y'); zlabel('Z');            
        end

        % Button pushed function: ClearFixationButton
        function ClearFixationButtonPushed(app, event)
            global fixingCond_;
            global F_; F_ = [];
            fixingCond_ = [];
            app.StressSimulationonSolidDomainPanel.Enable = 'off';
            app.MaterialLayoutDesignPanel.Enable = 'off';
            app.AdditionalStiffnessEvaluationPanel.Enable = 'off';
            ShowProblemDescriptionMenuSelected(app, event);
        end

        % Value changed function: SimulationTasksDropDown
        function SimulationTasksDropDownValueChanged(app, event)
            global cartesianStressField_;
            value = app.SimulationTasksDropDown.Value;
            switch value
                case 'None'
                    return;
                case 'Mtd - Topology Optimization'
                    app.comp_SimTask_TopologyOptimization_Func = Mtd_TopologyOptimization(app);
                case 'Mtd - Porous Infill Optimization'
                    app.comp_SimTask_TopologyOptimization_Func = Mtd_PorousInfillOptimization(app);              
                case 'Mtd - PSLs-guided Infill Design'
                    if isempty(cartesianStressField_)
                        warning('None Stress Field Existing!');
                        app.SimulationTasksDropDown.Value = 'None'; return;
                    end
                    app.comp_SimTask_PSLsGuidedStructDesign_Func = Mtd_PSLsGuidedStructDesign(app);
                case 'Mtd - Stress-aligned Volumetric Michell Trusses Infill Design'
                    if isempty(cartesianStressField_)
                        warning('None Stress Field Existing!');
                        app.SimulationTasksDropDown.Value = 'None'; return;
                    end
                    app.comp_SimTask_StressAwareGeoSyn3rdParty_AroraMethod = Mtd_StressAlignedVolumetricMichellTrusses(app);                   
                case 'Mtd - Stress-aligned Conforming Lattice Infill Design'
                    if isempty(cartesianStressField_)
                        warning('None Stress Field Existing!');
                        app.SimulationTasksDropDown.Value = 'None'; return;
                    end
                    app.comp_SimTask_StressAwareGeoSyn3rdParty_GaoMethod = Mtd_StressAlignedConformingLattice(app);                   
                case 'Mtd - Stress-aware Graded Voronoi Diagram Infill Design'
                    if isempty(cartesianStressField_)
                        warning('None Stress Field Existing!');
                        app.SimulationTasksDropDown.Value = 'None'; return;
                    end
                    app.comp_SimTask_StressAwareGeoSyn3rdParty_VoronoiMethod = Mtd_StressAwareGradedVoronoi(app);                                    
            end
            MainWindowCtrl(app, 0);
            app.VisualizationMenu.Enable = 'on';
            app.FileMenu.Enable = 'on';
        end

        % Menu selected function: ShowVertexEdgeGraphMenu
        function InputEdgeVertexGraphMenuSelected(app, event)
            global axHandle_;
            global vertexEdgeGraph_;
            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            [az, el] = view(axHandle_);
            cla(axHandle_); colorbar(axHandle_, 'off');
            Vis_DrawGraph3D(axHandle_, vertexEdgeGraph_.nodeCoords, vertexEdgeGraph_.eNodMat);
            view(axHandle_, az, el);
            Vis_UserLighting(axHandle_);
        end

        % Menu selected function: ShowDeformationMenu
        function TotalMenu_2Selected(app, event)
            global axHandle_;
            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            [az, el] = view(axHandle_);
            cla(axHandle_);            
            Vis_ShowDeformation(axHandle_, 'Total');
            view(axHandle_, az, el);             
        end

        % Menu selected function: ShowStressFieldvonMisesStressMenu
        function VonMisesMenuSelected(app, event)
            global axHandle_;
            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            [az, el] = view(axHandle_);
            cla(axHandle_);            
            Vis_ShowScalarStressComponents(axHandle_, 'Von Mises', 0);
            view(axHandle_, az, el);            
        end

        % Menu selected function: ShowPSLsMenu
        function PSLsMenuSelected(app, event)
            global axHandle_;
            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            [az, el] = view(axHandle_);
            cla(axHandle_); colorbar(axHandle_, 'off');           
            Vis_ShowPSLs(axHandle_);
            view(axHandle_, az, el);
            Vis_UserLighting(axHandle_);
        end

        % Button pushed function: FEAwithSolidDesignDomainButton
        function FEAwithSolidDesignDomainButtonPushed(app, event)
            global outPath_;
            global cartesianStressField_;
	        global vonMisesStressField_;
            global densityLayout_;
            global complianceSolid_;

            MainWindowCtrl(app, 0);
            app.LinearSystemSolverPanel.Enable = 'off';
            pause(1);

            GatherLSSandMPsettings(app);
            [complianceSolid_, ~] = FEA_ComputeComplianceVoxel();

            disp('Stress Analysis on Solid Domain ...');
            [cartesianStressField_, vonMisesStressField_] = FEA_StressAnalysis();
            vonMisesStressPerElement = FEA_ComputePerElementVonMisesStress(cartesianStressField_);
            dominantDirSolid = Common_ExtractDominantDirectionsFromPrincipalStressDirections(cartesianStressField_);
            niftiwrite(dominantDirSolid, strcat(outPath_, 'dominantDirSolid.nii'));
            vonMisesVolume = Common_ConvertPerEleVector2Volume(vonMisesStressPerElement);
            densityLayout = densityLayout_;
            densityLayout_ = ones(size(dominantDirSolid,1),1);
            IO_ExportDesignWithOneProperty_nii(vonMisesVolume, strcat(outPath_, 'ResultVolume_Solid_vonMises.nii'));
            densityLayout_ = densityLayout; 

            MainWindowCtrl(app, 1);            
            app.CellSizeEditField.Editable = 'off';
            app.LinearSystemSolverPanel.Enable = 'on';           
            app.SolidComplianceEditField.Value = complianceSolid_;
            
            app.ShowDeformationMenu.Enable = 'on'; 
            app.ShowStressFieldvonMisesStressMenu.Enable = 'on';
            pause(1);
            VonMisesMenuSelected(app);

            app.StiffnesswrtChangedLoadingDirectionsButton.Enable = 'on';
        end

        % Menu selected function: ShowDesignbyIsosurfaceNotrecommendedMenu
        function ShowDesignbyIsosurfaceNotrecommendedMenuSelected(app, event)
            global axHandle_;
            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            [az, el] = view(axHandle_);
            cla(axHandle_); colorbar(axHandle_, 'off');           
            Vis_ShowDesignByDensityLayoutInIsosurface(axHandle_);
            view(axHandle_, az, el);
            Vis_UserLighting(axHandle_);
        end

        % Menu selected function: ShowDesignDomainMenu
        function ShowDesignDomainMenuSelected(app, event)
            global axHandle_;
            if ~isvalid(axHandle_), axHandle_ = gca; view(axHandle_,3); end
            [az, el] = view(axHandle_);
            cla(axHandle_); colorbar(axHandle_, 'off');           
            Vis_ShowDesignDomain(axHandle_)
            view(axHandle_, az, el);
            Vis_UserLighting(axHandle_);            
        end

        % Menu selected function: VoxelModelTopVoxelMenu_2
        function ExportVoxelModelTopVoxelMenuSelected(app, event)
            [fileName, dataPath] = uiputfile('*.TopVoxel', 'Select a Path to Write');
            if isnumeric(fileName) || isnumeric(dataPath), return; end
            ofileName = strcat(dataPath,fileName);
            IO_ExportTopVoxels(ofileName);
        end

        % Menu selected function: DesignVolumeniiMenu
        function ExportDesigninVolumeniiMenuSelected(app, event)
            [fileName, dataPath] = uiputfile('*.nii', 'Select a Path to Write');
            if isnumeric(fileName) || isnumeric(dataPath), return; end
            ofileName = strcat(dataPath,fileName);
            IO_ExportDesignInVolume_nii(ofileName);               
        end

        % Value changed function: CenterXEditField, CenterYEditField, 
        % ...and 9 other components
        function SelectionOptionValueChanged(app, event)
            global hdSelectionBox_;
            value = app.SelectionOptionsDropDown.Value;            
            switch value
                case 'None'
                    if isempty(hdSelectionBox_), return; end
                    if isvalid(hdSelectionBox_)
                        set(hdSelectionBox_, 'visible', 'off');
                    end
                case 'Box'
                    cP1 = zeros(1,3); cP2 = cP1;
                    cP1(1) = app.CornerPot1XEditField.Value; cP1(2) = app.CornerPot1YEditField.Value; cP1(3) = app.CornerPot1ZEditField.Value;
                    cP2(1) = app.CornerPot2XEditField.Value; cP2(2) = app.CornerPot2YEditField.Value; cP2(3) = app.CornerPot2ZEditField.Value;
                    Vis_ShowSelectionBox(cP1, cP2);
                case 'Sphere'
                    sphereRad = abs(app.RadiusEditField.Value);
                    if 0==sphereRad, return; end                    
                    sphereCtr = zeros(1,3);
                    sphereCtr(1) = app.CenterXEditField.Value;
                    sphereCtr(2) = app.CenterYEditField.Value;
                    sphereCtr(3) = app.CenterZEditField.Value;
                    Vis_ShowSelectionSphere(sphereCtr, sphereRad);
            end
        end

        % Button pushed function: NodeSelectionButton_2
        function NodeSelectionButton_2Pushed(app, event)
            global axHandle_;
            SelectionOptionValueChanged(app, event);
            value = app.SelectionOptionsDropDown.Value;
            switch value
                case 'Box'
                    cP1 = zeros(1,3); cP2 = cP1;
                    cP1(1) = app.CornerPot1XEditField.Value; cP1(2) = app.CornerPot1YEditField.Value; cP1(3) = app.CornerPot1ZEditField.Value;
                    cP2(1) = app.CornerPot2XEditField.Value; cP2(2) = app.CornerPot2YEditField.Value; cP2(3) = app.CornerPot2ZEditField.Value;            
                    Interaction_PickBySelectionBox(axHandle_, cP1, cP2);
                case 'Sphere'
                    sphereRad = abs(app.RadiusEditField.Value);
                    if 0==sphereRad, return; end                    
                    sphereCtr = zeros(1,3);
                    sphereCtr(1) = app.CenterXEditField.Value;
                    sphereCtr(2) = app.CenterYEditField.Value;
                    sphereCtr(3) = app.CenterZEditField.Value;
                    Interaction_PickBySelectionShpere(axHandle_, sphereCtr, sphereRad);                    
            end
        end

        % Button pushed function: NodeUnSelectionButton
        function NodeUnSelectionButtonPushed(app, event)
            global axHandle_;
            SelectionOptionValueChanged(app, event);
            value = app.SelectionOptionsDropDown.Value;
            switch value
                case 'Box'
                    cP1 = zeros(1,3); cP2 = cP1;
                    cP1(1) = app.CornerPot1XEditField.Value; cP1(2) = app.CornerPot1YEditField.Value; cP1(3) = app.CornerPot1ZEditField.Value;
                    cP2(1) = app.CornerPot2XEditField.Value; cP2(2) = app.CornerPot2YEditField.Value; cP2(3) = app.CornerPot2ZEditField.Value;            
                    Interaction_UnPickBySelectionBox(axHandle_, cP1, cP2);   
                case 'Sphere'
                    sphereRad = abs(app.RadiusEditField.Value);
                    if 0==sphereRad, return; end
                    sphereCtr = zeros(1,3);
                    sphereCtr(1) = app.CenterXEditField.Value;
                    sphereCtr(2) = app.CenterYEditField.Value;
                    sphereCtr(3) = app.CenterZEditField.Value;
                    Interaction_UnPickBySelectionShpere(axHandle_, sphereCtr, sphereRad);  
            end         
        end

        % Menu selected function: ShowDesignVolumeWebGLMenu
        function ShowDesignVolumeWebGLMenuSelected(app, event)
            web('https://keksboter.github.io/quokka/');
        end

        % Menu selected function: VoxelModelwithDensityLayoutTopVoxelMenu
        function VoxelModelwithDensityLayoutTopVoxelMenuSelected(app, event)
            [fileName, dataPath] = uiputfile('*.TopVoxel', 'Select a Path to Write');
            if isnumeric(fileName) || isnumeric(dataPath), return; end
            ofileName = strcat(dataPath,fileName);
            IO_ExportTopVoxels(ofileName,1);            
        end

        % Button pushed function: 
        % StiffnessofExternalMeshGraphviaVoxelFEAButton
        function StiffnessofExternalMeshGraphviaVoxelFEAButtonPushed(app, event)
            app.comp_SimTask_MeshGraphBasedStructDesign_Func = Mtd_MeshGraphBasedStructDesign(app);
            MainWindowCtrl(app, 0);
            app.VisualizationMenu.Enable = 'on';
            app.FileMenu.Enable = 'on';            
        end

        % Button pushed function: 
        % StiffnesswrtChangedLoadingDirectionsButton
        function StiffnesswrtChangedLoadingDirectionsButtonPushed(app, event)
            global densityLayout_;
            if isempty(densityLayout_)
                warning('No Design is Available! The Solid Domain will be Taken!');
            end
            MainWindowCtrl(app, 0);
            ShowProblemDescriptionMenuSelected(app);
            app.comp_SimTask_DesignWRTchangedLoadingDirections_Func = Mtd_StructuralRobustnessExploration(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 404 768];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create FileMenu
            app.FileMenu = uimenu(app.UIFigure);
            app.FileMenu.Text = 'File';

            % Create ImportMenu
            app.ImportMenu = uimenu(app.FileMenu);
            app.ImportMenu.Text = 'Import';

            % Create TriangularSurfaceMeshplyobjstlMenu
            app.TriangularSurfaceMeshplyobjstlMenu = uimenu(app.ImportMenu);
            app.TriangularSurfaceMeshplyobjstlMenu.MenuSelectedFcn = createCallbackFcn(app, @SurfaceTriMeshplyobjMenuSelected, true);
            app.TriangularSurfaceMeshplyobjstlMenu.Text = 'Triangular Surface Mesh (*.ply, *.obj, *.stl)';

            % Create VoxelModelTopVoxelMenu
            app.VoxelModelTopVoxelMenu = uimenu(app.ImportMenu);
            app.VoxelModelTopVoxelMenu.MenuSelectedFcn = createCallbackFcn(app, @ImportVoxelModelTopVoxelMenuSelected, true);
            app.VoxelModelTopVoxelMenu.Text = 'Voxel Model (*.TopVoxel)';

            % Create BuiltinShapesMenu
            app.BuiltinShapesMenu = uimenu(app.ImportMenu);
            app.BuiltinShapesMenu.Text = 'Built-in Shapes';

            % Create CuboidMenu
            app.CuboidMenu = uimenu(app.BuiltinShapesMenu);
            app.CuboidMenu.MenuSelectedFcn = createCallbackFcn(app, @CuboidMenuSelected, true);
            app.CuboidMenu.Text = 'Cuboid';

            % Create LshapeMenu
            app.LshapeMenu = uimenu(app.BuiltinShapesMenu);
            app.LshapeMenu.MenuSelectedFcn = createCallbackFcn(app, @LshapeMenuSelected, true);
            app.LshapeMenu.Text = 'L-shape';

            % Create CylinderMenu
            app.CylinderMenu = uimenu(app.BuiltinShapesMenu);
            app.CylinderMenu.MenuSelectedFcn = createCallbackFcn(app, @CylinderMenuSelected, true);
            app.CylinderMenu.Text = 'Cylinder';

            % Create ExportMenu
            app.ExportMenu = uimenu(app.FileMenu);
            app.ExportMenu.Text = 'Export';

            % Create DesignVolumeniiMenu
            app.DesignVolumeniiMenu = uimenu(app.ExportMenu);
            app.DesignVolumeniiMenu.MenuSelectedFcn = createCallbackFcn(app, @ExportDesigninVolumeniiMenuSelected, true);
            app.DesignVolumeniiMenu.Text = 'Design Volume (*.nii)';

            % Create VoxelModelTopVoxelMenu_2
            app.VoxelModelTopVoxelMenu_2 = uimenu(app.ExportMenu);
            app.VoxelModelTopVoxelMenu_2.MenuSelectedFcn = createCallbackFcn(app, @ExportVoxelModelTopVoxelMenuSelected, true);
            app.VoxelModelTopVoxelMenu_2.Text = 'Voxel Model (*.TopVoxel)';

            % Create VoxelModelwithDensityLayoutTopVoxelMenu
            app.VoxelModelwithDensityLayoutTopVoxelMenu = uimenu(app.ExportMenu);
            app.VoxelModelwithDensityLayoutTopVoxelMenu.MenuSelectedFcn = createCallbackFcn(app, @VoxelModelwithDensityLayoutTopVoxelMenuSelected, true);
            app.VoxelModelwithDensityLayoutTopVoxelMenu.Text = 'Voxel Model with Density Layout (*.TopVoxel)';

            % Create VisualizationMenu
            app.VisualizationMenu = uimenu(app.UIFigure);
            app.VisualizationMenu.ForegroundColor = [0.149 0.149 0.149];
            app.VisualizationMenu.Text = 'Visualization';

            % Create ShowDesignVolumeWebGLMenu
            app.ShowDesignVolumeWebGLMenu = uimenu(app.VisualizationMenu);
            app.ShowDesignVolumeWebGLMenu.MenuSelectedFcn = createCallbackFcn(app, @ShowDesignVolumeWebGLMenuSelected, true);
            app.ShowDesignVolumeWebGLMenu.Text = 'Show Design Volume (WebGL)';

            % Create ShowInputTriangularSurfaceMeshMenu
            app.ShowInputTriangularSurfaceMeshMenu = uimenu(app.VisualizationMenu);
            app.ShowInputTriangularSurfaceMeshMenu.MenuSelectedFcn = createCallbackFcn(app, @ShowInputTriSurfaceMeshMenuSelected, true);
            app.ShowInputTriangularSurfaceMeshMenu.Text = 'Show Input Triangular Surface Mesh';

            % Create ShowProblemDescriptionMenu
            app.ShowProblemDescriptionMenu = uimenu(app.VisualizationMenu);
            app.ShowProblemDescriptionMenu.MenuSelectedFcn = createCallbackFcn(app, @ShowProblemDescriptionMenuSelected, true);
            app.ShowProblemDescriptionMenu.Text = 'Show Problem Description';

            % Create ShowDesignDomainMenu
            app.ShowDesignDomainMenu = uimenu(app.VisualizationMenu);
            app.ShowDesignDomainMenu.MenuSelectedFcn = createCallbackFcn(app, @ShowDesignDomainMenuSelected, true);
            app.ShowDesignDomainMenu.Text = 'Show Design Domain';

            % Create ShowDeformationMenu
            app.ShowDeformationMenu = uimenu(app.VisualizationMenu);
            app.ShowDeformationMenu.MenuSelectedFcn = createCallbackFcn(app, @TotalMenu_2Selected, true);
            app.ShowDeformationMenu.Text = 'Show Deformation';

            % Create ShowStressFieldvonMisesStressMenu
            app.ShowStressFieldvonMisesStressMenu = uimenu(app.VisualizationMenu);
            app.ShowStressFieldvonMisesStressMenu.MenuSelectedFcn = createCallbackFcn(app, @VonMisesMenuSelected, true);
            app.ShowStressFieldvonMisesStressMenu.Text = 'Show Stress Field (von Mises Stress)';

            % Create ShowPSLsMenu
            app.ShowPSLsMenu = uimenu(app.VisualizationMenu);
            app.ShowPSLsMenu.MenuSelectedFcn = createCallbackFcn(app, @PSLsMenuSelected, true);
            app.ShowPSLsMenu.Text = 'Show PSLs';

            % Create ShowVertexEdgeGraphMenu
            app.ShowVertexEdgeGraphMenu = uimenu(app.VisualizationMenu);
            app.ShowVertexEdgeGraphMenu.MenuSelectedFcn = createCallbackFcn(app, @InputEdgeVertexGraphMenuSelected, true);
            app.ShowVertexEdgeGraphMenu.Text = 'Show Vertex-Edge Graph';

            % Create ShowDesignbyIsosurfaceNotrecommendedMenu
            app.ShowDesignbyIsosurfaceNotrecommendedMenu = uimenu(app.VisualizationMenu);
            app.ShowDesignbyIsosurfaceNotrecommendedMenu.MenuSelectedFcn = createCallbackFcn(app, @ShowDesignbyIsosurfaceNotrecommendedMenuSelected, true);
            app.ShowDesignbyIsosurfaceNotrecommendedMenu.Text = 'Show Design by Isosurface (Not recommended)';

            % Create TabGroup3
            app.TabGroup3 = uitabgroup(app.UIFigure);
            app.TabGroup3.Position = [0 5 403 762];

            % Create ModelingTab
            app.ModelingTab = uitab(app.TabGroup3);
            app.ModelingTab.Title = 'Modeling';

            % Create DomainVoxelizationPanel
            app.DomainVoxelizationPanel = uipanel(app.ModelingTab);
            app.DomainVoxelizationPanel.Title = 'Domain Voxelization';
            app.DomainVoxelizationPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.DomainVoxelizationPanel.FontWeight = 'bold';
            app.DomainVoxelizationPanel.Position = [0 506 400 232];

            % Create TargetVoxelResolutionEditFieldLabel
            app.TargetVoxelResolutionEditFieldLabel = uilabel(app.DomainVoxelizationPanel);
            app.TargetVoxelResolutionEditFieldLabel.HorizontalAlignment = 'right';
            app.TargetVoxelResolutionEditFieldLabel.Position = [135 180 131 22];
            app.TargetVoxelResolutionEditFieldLabel.Text = 'Target Voxel Resolution';

            % Create TargetVoxelResolutionEditField
            app.TargetVoxelResolutionEditField = uieditfield(app.DomainVoxelizationPanel, 'numeric');
            app.TargetVoxelResolutionEditField.ValueDisplayFormat = '%.0f';
            app.TargetVoxelResolutionEditField.Position = [281 180 100 22];
            app.TargetVoxelResolutionEditField.Value = 256;

            % Create CellSizeEditFieldLabel
            app.CellSizeEditFieldLabel = uilabel(app.DomainVoxelizationPanel);
            app.CellSizeEditFieldLabel.HorizontalAlignment = 'right';
            app.CellSizeEditFieldLabel.Position = [214 139 52 22];
            app.CellSizeEditFieldLabel.Text = 'Cell Size';

            % Create CellSizeEditField
            app.CellSizeEditField = uieditfield(app.DomainVoxelizationPanel, 'numeric');
            app.CellSizeEditField.Position = [281 139 100 22];
            app.CellSizeEditField.Value = 1;

            % Create VoxelizingButton
            app.VoxelizingButton = uibutton(app.DomainVoxelizationPanel, 'push');
            app.VoxelizingButton.ButtonPushedFcn = createCallbackFcn(app, @VoxelizingButtonPushed, true);
            app.VoxelizingButton.FontWeight = 'bold';
            app.VoxelizingButton.Position = [282 97 100 23];
            app.VoxelizingButton.Text = 'Voxelizing';

            % Create ElementsEditFieldLabel
            app.ElementsEditFieldLabel = uilabel(app.DomainVoxelizationPanel);
            app.ElementsEditFieldLabel.HorizontalAlignment = 'right';
            app.ElementsEditFieldLabel.Position = [209 54 62 22];
            app.ElementsEditFieldLabel.Text = '#Elements';

            % Create ElementsEditField
            app.ElementsEditField = uieditfield(app.DomainVoxelizationPanel, 'numeric');
            app.ElementsEditField.ValueDisplayFormat = '%.0f';
            app.ElementsEditField.Position = [286 54 96 22];

            % Create DOFsEditFieldLabel
            app.DOFsEditFieldLabel = uilabel(app.DomainVoxelizationPanel);
            app.DOFsEditFieldLabel.HorizontalAlignment = 'right';
            app.DOFsEditFieldLabel.Position = [228 12 43 22];
            app.DOFsEditFieldLabel.Text = '#DOFs';

            % Create DOFsEditField
            app.DOFsEditField = uieditfield(app.DomainVoxelizationPanel, 'numeric');
            app.DOFsEditField.ValueDisplayFormat = '%.0f';
            app.DOFsEditField.Position = [286 12 97 22];

            % Create ApplyforBoundaryConditionsPanel
            app.ApplyforBoundaryConditionsPanel = uipanel(app.ModelingTab);
            app.ApplyforBoundaryConditionsPanel.Title = 'Apply for Boundary Conditions';
            app.ApplyforBoundaryConditionsPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.ApplyforBoundaryConditionsPanel.FontWeight = 'bold';
            app.ApplyforBoundaryConditionsPanel.Position = [0 -6 400 512];

            % Create NodeUnSelectionButton
            app.NodeUnSelectionButton = uibutton(app.ApplyforBoundaryConditionsPanel, 'push');
            app.NodeUnSelectionButton.ButtonPushedFcn = createCallbackFcn(app, @NodeUnSelectionButtonPushed, true);
            app.NodeUnSelectionButton.BackgroundColor = [0.9608 0.9608 0.9608];
            app.NodeUnSelectionButton.FontWeight = 'bold';
            app.NodeUnSelectionButton.Position = [42 276 122 23];
            app.NodeUnSelectionButton.Text = 'Node Un-Selection';

            % Create NodeSelectionButton_2
            app.NodeSelectionButton_2 = uibutton(app.ApplyforBoundaryConditionsPanel, 'push');
            app.NodeSelectionButton_2.ButtonPushedFcn = createCallbackFcn(app, @NodeSelectionButton_2Pushed, true);
            app.NodeSelectionButton_2.BackgroundColor = [0.9608 0.9608 0.9608];
            app.NodeSelectionButton_2.FontWeight = 'bold';
            app.NodeSelectionButton_2.Position = [242 276 122 23];
            app.NodeSelectionButton_2.Text = 'Node Selection';

            % Create TabGroup2
            app.TabGroup2 = uitabgroup(app.ApplyforBoundaryConditionsPanel);
            app.TabGroup2.Position = [1 84 400 176];

            % Create FixingTab
            app.FixingTab = uitab(app.TabGroup2);
            app.FixingTab.Title = 'Fixing';
            app.FixingTab.BackgroundColor = [0.9412 0.9412 0.9412];

            % Create ZDirFixedCheckBox
            app.ZDirFixedCheckBox = uicheckbox(app.FixingTab);
            app.ZDirFixedCheckBox.Text = 'Z-Dir Fixed';
            app.ZDirFixedCheckBox.Position = [311 49 81 22];
            app.ZDirFixedCheckBox.Value = true;

            % Create ApplyforFixationButton
            app.ApplyforFixationButton = uibutton(app.FixingTab, 'push');
            app.ApplyforFixationButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyforFixationButtonPushed, true);
            app.ApplyforFixationButton.FontWeight = 'bold';
            app.ApplyforFixationButton.Position = [227 13 165 26];
            app.ApplyforFixationButton.Text = 'Apply for Fixation';

            % Create ClearFixationButton
            app.ClearFixationButton = uibutton(app.FixingTab, 'push');
            app.ClearFixationButton.ButtonPushedFcn = createCallbackFcn(app, @ClearFixationButtonPushed, true);
            app.ClearFixationButton.FontWeight = 'bold';
            app.ClearFixationButton.Position = [89 15 100 23];
            app.ClearFixationButton.Text = 'Clear Fixation';

            % Create YDirFixedCheckBox
            app.YDirFixedCheckBox = uicheckbox(app.FixingTab);
            app.YDirFixedCheckBox.Text = 'Y-Dir Fixed';
            app.YDirFixedCheckBox.Position = [311 83 81 22];
            app.YDirFixedCheckBox.Value = true;

            % Create XDirFixedCheckBox
            app.XDirFixedCheckBox = uicheckbox(app.FixingTab);
            app.XDirFixedCheckBox.Text = 'X-Dir Fixed';
            app.XDirFixedCheckBox.Position = [311 120 82 22];
            app.XDirFixedCheckBox.Value = true;

            % Create LoadingTab
            app.LoadingTab = uitab(app.TabGroup2);
            app.LoadingTab.Title = 'Loading';
            app.LoadingTab.BackgroundColor = [0.9412 0.9412 0.9412];

            % Create FxLabel
            app.FxLabel = uilabel(app.LoadingTab);
            app.FxLabel.Position = [234 120 38 22];
            app.FxLabel.Text = 'Fx (N)';

            % Create FxNEditField
            app.FxNEditField = uieditfield(app.LoadingTab, 'numeric');
            app.FxNEditField.Position = [274 120 113 22];

            % Create ApplyforLoadsButton
            app.ApplyforLoadsButton = uibutton(app.LoadingTab, 'push');
            app.ApplyforLoadsButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyforLoadsButtonPushed, true);
            app.ApplyforLoadsButton.FontWeight = 'bold';
            app.ApplyforLoadsButton.Position = [238 13 150 26];
            app.ApplyforLoadsButton.Text = 'Apply for Loads';

            % Create ClearLoadsButton
            app.ClearLoadsButton = uibutton(app.LoadingTab, 'push');
            app.ClearLoadsButton.ButtonPushedFcn = createCallbackFcn(app, @ClearLoadsButtonPushed, true);
            app.ClearLoadsButton.FontWeight = 'bold';
            app.ClearLoadsButton.Position = [87 15 100 23];
            app.ClearLoadsButton.Text = 'Clear Loads';

            % Create FyNEditFieldLabel
            app.FyNEditFieldLabel = uilabel(app.LoadingTab);
            app.FyNEditFieldLabel.Position = [235 83 38 22];
            app.FyNEditFieldLabel.Text = 'Fy (N)';

            % Create FyNEditField
            app.FyNEditField = uieditfield(app.LoadingTab, 'numeric');
            app.FyNEditField.Position = [275 83 113 22];

            % Create FzNEditFieldLabel
            app.FzNEditFieldLabel = uilabel(app.LoadingTab);
            app.FzNEditFieldLabel.Position = [235 49 38 22];
            app.FzNEditFieldLabel.Text = 'Fz (N)';

            % Create FzNEditField
            app.FzNEditField = uieditfield(app.LoadingTab, 'numeric');
            app.FzNEditField.Position = [275 49 113 22];

            % Create BuiltinBoundaryConditionsPanel
            app.BuiltinBoundaryConditionsPanel = uipanel(app.ApplyforBoundaryConditionsPanel);
            app.BuiltinBoundaryConditionsPanel.Title = 'Built-in Boundary Conditions';
            app.BuiltinBoundaryConditionsPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.BuiltinBoundaryConditionsPanel.Position = [1 1 400 84];

            % Create OnlyCuboidDomainDropDownLabel
            app.OnlyCuboidDomainDropDownLabel = uilabel(app.BuiltinBoundaryConditionsPanel);
            app.OnlyCuboidDomainDropDownLabel.HorizontalAlignment = 'right';
            app.OnlyCuboidDomainDropDownLabel.Position = [161 20 116 22];
            app.OnlyCuboidDomainDropDownLabel.Text = 'Only Cuboid Domain';

            % Create OnlyCuboidDomainDropDown
            app.OnlyCuboidDomainDropDown = uidropdown(app.BuiltinBoundaryConditionsPanel);
            app.OnlyCuboidDomainDropDown.Items = {'None', 'Cuboid - Cantilever 1', 'Cuboid - Cantilever 2', 'Cuboid - Cantilever 3', 'Cuboid - Cantilever 4', 'Cuboid - MBB Half', 'Cuboid - MBB'};
            app.OnlyCuboidDomainDropDown.ValueChangedFcn = createCallbackFcn(app, @OnlyCuboidDomainDropDownValueChanged, true);
            app.OnlyCuboidDomainDropDown.Position = [284 20 100 22];
            app.OnlyCuboidDomainDropDown.Value = 'None';

            % Create SelectionOptionsDropDownLabel
            app.SelectionOptionsDropDownLabel = uilabel(app.ApplyforBoundaryConditionsPanel);
            app.SelectionOptionsDropDownLabel.HorizontalAlignment = 'right';
            app.SelectionOptionsDropDownLabel.Position = [171 460 99 22];
            app.SelectionOptionsDropDownLabel.Text = 'Selection Options';

            % Create SelectionOptionsDropDown
            app.SelectionOptionsDropDown = uidropdown(app.ApplyforBoundaryConditionsPanel);
            app.SelectionOptionsDropDown.Items = {'None', 'Box', 'Sphere'};
            app.SelectionOptionsDropDown.ValueChangedFcn = createCallbackFcn(app, @SelectionOptionValueChanged, true);
            app.SelectionOptionsDropDown.Position = [285 460 100 22];
            app.SelectionOptionsDropDown.Value = 'Box';

            % Create TabGroupSelection
            app.TabGroupSelection = uitabgroup(app.ApplyforBoundaryConditionsPanel);
            app.TabGroupSelection.Position = [0 310 399 142];

            % Create BoxSelectionTab
            app.BoxSelectionTab = uitab(app.TabGroupSelection);
            app.BoxSelectionTab.Title = 'Box Selection';

            % Create CornerPot1XLabel
            app.CornerPot1XLabel = uilabel(app.BoxSelectionTab);
            app.CornerPot1XLabel.HorizontalAlignment = 'right';
            app.CornerPot1XLabel.Position = [17 70 104 22];
            app.CornerPot1XLabel.Text = 'Corner Pot 1 (*): X';

            % Create CornerPot1XEditField
            app.CornerPot1XEditField = uieditfield(app.BoxSelectionTab, 'numeric');
            app.CornerPot1XEditField.ValueChangedFcn = createCallbackFcn(app, @SelectionOptionValueChanged, true);
            app.CornerPot1XEditField.Position = [136 70 39 22];

            % Create CornerPot2XLabel
            app.CornerPot2XLabel = uilabel(app.BoxSelectionTab);
            app.CornerPot2XLabel.HorizontalAlignment = 'right';
            app.CornerPot2XLabel.Position = [222 70 106 22];
            app.CornerPot2XLabel.Text = 'Corner Pot 2 (+): X';

            % Create CornerPot2XEditField
            app.CornerPot2XEditField = uieditfield(app.BoxSelectionTab, 'numeric');
            app.CornerPot2XEditField.ValueChangedFcn = createCallbackFcn(app, @SelectionOptionValueChanged, true);
            app.CornerPot2XEditField.Position = [343 70 39 22];

            % Create CornerPot1YEditFieldLabel
            app.CornerPot1YEditFieldLabel = uilabel(app.BoxSelectionTab);
            app.CornerPot1YEditFieldLabel.HorizontalAlignment = 'right';
            app.CornerPot1YEditFieldLabel.Position = [19 39 103 22];
            app.CornerPot1YEditFieldLabel.Text = 'Corner Pot 1 (*): Y';

            % Create CornerPot1YEditField
            app.CornerPot1YEditField = uieditfield(app.BoxSelectionTab, 'numeric');
            app.CornerPot1YEditField.ValueChangedFcn = createCallbackFcn(app, @SelectionOptionValueChanged, true);
            app.CornerPot1YEditField.Position = [137 39 39 22];

            % Create CornerPot2YLabel
            app.CornerPot2YLabel = uilabel(app.BoxSelectionTab);
            app.CornerPot2YLabel.HorizontalAlignment = 'right';
            app.CornerPot2YLabel.Position = [223 39 106 22];
            app.CornerPot2YLabel.Text = 'Corner Pot 2 (+): Y';

            % Create CornerPot2YEditField
            app.CornerPot2YEditField = uieditfield(app.BoxSelectionTab, 'numeric');
            app.CornerPot2YEditField.ValueChangedFcn = createCallbackFcn(app, @SelectionOptionValueChanged, true);
            app.CornerPot2YEditField.Position = [343 39 39 22];

            % Create CornerPot1ZEditFieldLabel
            app.CornerPot1ZEditFieldLabel = uilabel(app.BoxSelectionTab);
            app.CornerPot1ZEditFieldLabel.HorizontalAlignment = 'right';
            app.CornerPot1ZEditFieldLabel.Position = [19 8 103 22];
            app.CornerPot1ZEditFieldLabel.Text = 'Corner Pot 1 (*): Z';

            % Create CornerPot1ZEditField
            app.CornerPot1ZEditField = uieditfield(app.BoxSelectionTab, 'numeric');
            app.CornerPot1ZEditField.ValueChangedFcn = createCallbackFcn(app, @SelectionOptionValueChanged, true);
            app.CornerPot1ZEditField.Position = [137 8 39 22];

            % Create CornerPot2ZEditFieldLabel
            app.CornerPot2ZEditFieldLabel = uilabel(app.BoxSelectionTab);
            app.CornerPot2ZEditFieldLabel.HorizontalAlignment = 'right';
            app.CornerPot2ZEditFieldLabel.Position = [224 8 105 22];
            app.CornerPot2ZEditFieldLabel.Text = 'Corner Pot 2 (+): Z';

            % Create CornerPot2ZEditField
            app.CornerPot2ZEditField = uieditfield(app.BoxSelectionTab, 'numeric');
            app.CornerPot2ZEditField.ValueChangedFcn = createCallbackFcn(app, @SelectionOptionValueChanged, true);
            app.CornerPot2ZEditField.Position = [343 8 39 22];

            % Create SphereSelectionTab
            app.SphereSelectionTab = uitab(app.TabGroupSelection);
            app.SphereSelectionTab.Title = 'Sphere Selection';

            % Create CenterXEditFieldLabel
            app.CenterXEditFieldLabel = uilabel(app.SphereSelectionTab);
            app.CenterXEditFieldLabel.HorizontalAlignment = 'right';
            app.CenterXEditFieldLabel.Position = [13 78 52 22];
            app.CenterXEditFieldLabel.Text = 'Center X';

            % Create CenterXEditField
            app.CenterXEditField = uieditfield(app.SphereSelectionTab, 'numeric');
            app.CenterXEditField.ValueChangedFcn = createCallbackFcn(app, @SelectionOptionValueChanged, true);
            app.CenterXEditField.Position = [80 78 100 22];

            % Create CenterYEditFieldLabel
            app.CenterYEditFieldLabel = uilabel(app.SphereSelectionTab);
            app.CenterYEditFieldLabel.HorizontalAlignment = 'right';
            app.CenterYEditFieldLabel.Position = [13 49 52 22];
            app.CenterYEditFieldLabel.Text = 'Center Y';

            % Create CenterYEditField
            app.CenterYEditField = uieditfield(app.SphereSelectionTab, 'numeric');
            app.CenterYEditField.ValueChangedFcn = createCallbackFcn(app, @SelectionOptionValueChanged, true);
            app.CenterYEditField.Position = [80 49 100 22];

            % Create CenterZEditFieldLabel
            app.CenterZEditFieldLabel = uilabel(app.SphereSelectionTab);
            app.CenterZEditFieldLabel.HorizontalAlignment = 'right';
            app.CenterZEditFieldLabel.Position = [14 19 52 22];
            app.CenterZEditFieldLabel.Text = 'Center Z';

            % Create CenterZEditField
            app.CenterZEditField = uieditfield(app.SphereSelectionTab, 'numeric');
            app.CenterZEditField.ValueChangedFcn = createCallbackFcn(app, @SelectionOptionValueChanged, true);
            app.CenterZEditField.Position = [81 19 100 22];

            % Create RadiusEditFieldLabel
            app.RadiusEditFieldLabel = uilabel(app.SphereSelectionTab);
            app.RadiusEditFieldLabel.HorizontalAlignment = 'right';
            app.RadiusEditFieldLabel.Position = [225 48 42 22];
            app.RadiusEditFieldLabel.Text = 'Radius';

            % Create RadiusEditField
            app.RadiusEditField = uieditfield(app.SphereSelectionTab, 'numeric');
            app.RadiusEditField.ValueChangedFcn = createCallbackFcn(app, @SelectionOptionValueChanged, true);
            app.RadiusEditField.Position = [282 48 100 22];
            app.RadiusEditField.Value = 10;

            % Create SimulationTab
            app.SimulationTab = uitab(app.TabGroup3);
            app.SimulationTab.Title = 'Simulation';

            % Create MaterialPropertiesPanel
            app.MaterialPropertiesPanel = uipanel(app.SimulationTab);
            app.MaterialPropertiesPanel.Title = 'Material Properties';
            app.MaterialPropertiesPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.MaterialPropertiesPanel.Position = [0 571 400 167];

            % Create YoungsModulusStiffmaterialEditFieldLabel
            app.YoungsModulusStiffmaterialEditFieldLabel = uilabel(app.MaterialPropertiesPanel);
            app.YoungsModulusStiffmaterialEditFieldLabel.HorizontalAlignment = 'right';
            app.YoungsModulusStiffmaterialEditFieldLabel.Position = [96 104 173 22];
            app.YoungsModulusStiffmaterialEditFieldLabel.Text = 'Young''s Modulus (Stiff material)';

            % Create YoungsModulusStiffmaterialEditField
            app.YoungsModulusStiffmaterialEditField = uieditfield(app.MaterialPropertiesPanel, 'numeric');
            app.YoungsModulusStiffmaterialEditField.Position = [284 104 100 22];
            app.YoungsModulusStiffmaterialEditField.Value = 1;

            % Create PoissonsRatioEditFieldLabel
            app.PoissonsRatioEditFieldLabel = uilabel(app.MaterialPropertiesPanel);
            app.PoissonsRatioEditFieldLabel.HorizontalAlignment = 'right';
            app.PoissonsRatioEditFieldLabel.Position = [182 63 87 22];
            app.PoissonsRatioEditFieldLabel.Text = 'Poisson''s Ratio';

            % Create PoissonsRatioEditField
            app.PoissonsRatioEditField = uieditfield(app.MaterialPropertiesPanel, 'numeric');
            app.PoissonsRatioEditField.Position = [284 63 100 22];
            app.PoissonsRatioEditField.Value = 0.3;

            % Create YoungsModulusCompliantmaterialEditFieldLabel
            app.YoungsModulusCompliantmaterialEditFieldLabel = uilabel(app.MaterialPropertiesPanel);
            app.YoungsModulusCompliantmaterialEditFieldLabel.HorizontalAlignment = 'right';
            app.YoungsModulusCompliantmaterialEditFieldLabel.Position = [61 23 208 22];
            app.YoungsModulusCompliantmaterialEditFieldLabel.Text = 'Young''s Modulus (Compliant material)';

            % Create YoungsModulusCompliantmaterialEditField
            app.YoungsModulusCompliantmaterialEditField = uieditfield(app.MaterialPropertiesPanel, 'numeric');
            app.YoungsModulusCompliantmaterialEditField.Position = [284 23 100 22];
            app.YoungsModulusCompliantmaterialEditField.Value = 1e-06;

            % Create LinearSystemSolverPanel
            app.LinearSystemSolverPanel = uipanel(app.SimulationTab);
            app.LinearSystemSolverPanel.Title = 'Linear System Solver';
            app.LinearSystemSolverPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.LinearSystemSolverPanel.FontWeight = 'bold';
            app.LinearSystemSolverPanel.Position = [0 403 400 169];

            % Create ResidualThresholdEditFieldLabel
            app.ResidualThresholdEditFieldLabel = uilabel(app.LinearSystemSolverPanel);
            app.ResidualThresholdEditFieldLabel.HorizontalAlignment = 'right';
            app.ResidualThresholdEditFieldLabel.Position = [161 103 108 22];
            app.ResidualThresholdEditFieldLabel.Text = 'Residual Threshold';

            % Create ResidualThresholdEditField
            app.ResidualThresholdEditField = uieditfield(app.LinearSystemSolverPanel, 'numeric');
            app.ResidualThresholdEditField.Position = [284 103 100 22];
            app.ResidualThresholdEditField.Value = 0.001;

            % Create MaximumIterationsEditFieldLabel
            app.MaximumIterationsEditFieldLabel = uilabel(app.LinearSystemSolverPanel);
            app.MaximumIterationsEditFieldLabel.HorizontalAlignment = 'right';
            app.MaximumIterationsEditFieldLabel.Position = [153 62 116 22];
            app.MaximumIterationsEditFieldLabel.Text = '#Maximum Iterations';

            % Create MaximumIterationsEditField
            app.MaximumIterationsEditField = uieditfield(app.LinearSystemSolverPanel, 'numeric');
            app.MaximumIterationsEditField.ValueDisplayFormat = '%.0f';
            app.MaximumIterationsEditField.Position = [284 62 100 22];
            app.MaximumIterationsEditField.Value = 500;

            % Create WeightingFactorofJacobiSmoothingProcessEditFieldLabel
            app.WeightingFactorofJacobiSmoothingProcessEditFieldLabel = uilabel(app.LinearSystemSolverPanel);
            app.WeightingFactorofJacobiSmoothingProcessEditFieldLabel.HorizontalAlignment = 'right';
            app.WeightingFactorofJacobiSmoothingProcessEditFieldLabel.Position = [13 21 256 22];
            app.WeightingFactorofJacobiSmoothingProcessEditFieldLabel.Text = 'Weighting Factor of Jacobi Smoothing Process';

            % Create WeightingFactorofJacobiSmoothingProcessEditField
            app.WeightingFactorofJacobiSmoothingProcessEditField = uieditfield(app.LinearSystemSolverPanel, 'numeric');
            app.WeightingFactorofJacobiSmoothingProcessEditField.Position = [284 21 100 22];
            app.WeightingFactorofJacobiSmoothingProcessEditField.Value = 0.35;

            % Create MaterialLayoutDesignPanel
            app.MaterialLayoutDesignPanel = uipanel(app.SimulationTab);
            app.MaterialLayoutDesignPanel.Title = 'Material Layout Design';
            app.MaterialLayoutDesignPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.MaterialLayoutDesignPanel.FontWeight = 'bold';
            app.MaterialLayoutDesignPanel.Position = [0 127 400 101];

            % Create SimulationTasksDropDownLabel
            app.SimulationTasksDropDownLabel = uilabel(app.MaterialLayoutDesignPanel);
            app.SimulationTasksDropDownLabel.HorizontalAlignment = 'right';
            app.SimulationTasksDropDownLabel.FontWeight = 'bold';
            app.SimulationTasksDropDownLabel.Position = [82 32 103 22];
            app.SimulationTasksDropDownLabel.Text = 'Simulation Tasks';

            % Create SimulationTasksDropDown
            app.SimulationTasksDropDown = uidropdown(app.MaterialLayoutDesignPanel);
            app.SimulationTasksDropDown.Items = {'None', 'Mtd - Topology Optimization', 'Mtd - Porous Infill Optimization', 'Mtd - Stress-aware Graded Voronoi Diagram Infill Design', 'Mtd - PSLs-guided Infill Design', 'Mtd - Stress-aligned Conforming Lattice Infill Design', 'Mtd - Stress-aligned Volumetric Michell Trusses Infill Design'};
            app.SimulationTasksDropDown.ValueChangedFcn = createCallbackFcn(app, @SimulationTasksDropDownValueChanged, true);
            app.SimulationTasksDropDown.Position = [195 32 189 22];
            app.SimulationTasksDropDown.Value = 'None';

            % Create StressSimulationonSolidDomainPanel
            app.StressSimulationonSolidDomainPanel = uipanel(app.SimulationTab);
            app.StressSimulationonSolidDomainPanel.Title = 'Stress Simulation on Solid Domain';
            app.StressSimulationonSolidDomainPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.StressSimulationonSolidDomainPanel.FontWeight = 'bold';
            app.StressSimulationonSolidDomainPanel.Position = [0 229 400 174];

            % Create SolidComplianceEditFieldLabel
            app.SolidComplianceEditFieldLabel = uilabel(app.StressSimulationonSolidDomainPanel);
            app.SolidComplianceEditFieldLabel.HorizontalAlignment = 'right';
            app.SolidComplianceEditFieldLabel.Position = [211 70 98 22];
            app.SolidComplianceEditFieldLabel.Text = 'Solid Compliance';

            % Create SolidComplianceEditField
            app.SolidComplianceEditField = uieditfield(app.StressSimulationonSolidDomainPanel, 'numeric');
            app.SolidComplianceEditField.Position = [324 70 69 22];

            % Create FEAwithSolidDesignDomainButton
            app.FEAwithSolidDesignDomainButton = uibutton(app.StressSimulationonSolidDomainPanel, 'push');
            app.FEAwithSolidDesignDomainButton.ButtonPushedFcn = createCallbackFcn(app, @FEAwithSolidDesignDomainButtonPushed, true);
            app.FEAwithSolidDesignDomainButton.FontWeight = 'bold';
            app.FEAwithSolidDesignDomainButton.Position = [204 23 190 23];
            app.FEAwithSolidDesignDomainButton.Text = 'FEA with Solid Design Domain';

            % Create DesignComplianceEditFieldLabel
            app.DesignComplianceEditFieldLabel = uilabel(app.StressSimulationonSolidDomainPanel);
            app.DesignComplianceEditFieldLabel.HorizontalAlignment = 'right';
            app.DesignComplianceEditFieldLabel.Position = [6 70 109 22];
            app.DesignComplianceEditFieldLabel.Text = 'Design Compliance';

            % Create DesignComplianceEditField
            app.DesignComplianceEditField = uieditfield(app.StressSimulationonSolidDomainPanel, 'numeric');
            app.DesignComplianceEditField.Position = [130 70 69 22];

            % Create DesignVolEditFieldLabel
            app.DesignVolEditFieldLabel = uilabel(app.StressSimulationonSolidDomainPanel);
            app.DesignVolEditFieldLabel.HorizontalAlignment = 'right';
            app.DesignVolEditFieldLabel.Position = [47 119 66 22];
            app.DesignVolEditFieldLabel.Text = 'Design Vol.';

            % Create DesignVolEditField
            app.DesignVolEditField = uieditfield(app.StressSimulationonSolidDomainPanel, 'numeric');
            app.DesignVolEditField.Position = [128 119 69 22];

            % Create AdditionalStiffnessEvaluationPanel
            app.AdditionalStiffnessEvaluationPanel = uipanel(app.SimulationTab);
            app.AdditionalStiffnessEvaluationPanel.Title = 'Additional Stiffness Evaluation';
            app.AdditionalStiffnessEvaluationPanel.FontWeight = 'bold';
            app.AdditionalStiffnessEvaluationPanel.Position = [0 4 400 124];

            % Create StiffnessofExternalMeshGraphviaVoxelFEAButton
            app.StiffnessofExternalMeshGraphviaVoxelFEAButton = uibutton(app.AdditionalStiffnessEvaluationPanel, 'push');
            app.StiffnessofExternalMeshGraphviaVoxelFEAButton.ButtonPushedFcn = createCallbackFcn(app, @StiffnessofExternalMeshGraphviaVoxelFEAButtonPushed, true);
            app.StiffnessofExternalMeshGraphviaVoxelFEAButton.FontWeight = 'bold';
            app.StiffnessofExternalMeshGraphviaVoxelFEAButton.Position = [103 63 289 23];
            app.StiffnessofExternalMeshGraphviaVoxelFEAButton.Text = 'Stiffness of External Mesh/Graph via Voxel FEA';

            % Create StiffnesswrtChangedLoadingDirectionsButton
            app.StiffnesswrtChangedLoadingDirectionsButton = uibutton(app.AdditionalStiffnessEvaluationPanel, 'push');
            app.StiffnesswrtChangedLoadingDirectionsButton.ButtonPushedFcn = createCallbackFcn(app, @StiffnesswrtChangedLoadingDirectionsButtonPushed, true);
            app.StiffnesswrtChangedLoadingDirectionsButton.FontWeight = 'bold';
            app.StiffnesswrtChangedLoadingDirectionsButton.Position = [126 22 266 23];
            app.StiffnesswrtChangedLoadingDirectionsButton.Text = 'Stiffness w.r.t. Changed Loading Directions';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SGLDBench_Main(varargin)

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