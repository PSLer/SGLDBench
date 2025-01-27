classdef Mtd_TopologyOptimization < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        TopologyOptimizationPanel       matlab.ui.container.Panel
        ResultDisplayPanel              matlab.ui.container.Panel
        DesignVolumeFractionEditField   matlab.ui.control.NumericEditField
        DesignVolumeFractionEditFieldLabel  matlab.ui.control.Label
        DesignComplianceEditField       matlab.ui.control.NumericEditField
        DesignComplianceEditFieldLabel  matlab.ui.control.Label
        OptimizationProcessSettingsPanel  matlab.ui.container.Panel
        ConvergeMinChangeEditField      matlab.ui.control.NumericEditField
        ConvergeMinChangeEditFieldLabel  matlab.ui.control.Label
        ConvergeSharpnessEditField      matlab.ui.control.NumericEditField
        ConvergeSharpnessEditFieldLabel  matlab.ui.control.Label
        MovingStepsizeEditField         matlab.ui.control.NumericEditField
        MovingStepsizeEditField_2Label  matlab.ui.control.Label
        OptimizerDropDown               matlab.ui.control.DropDown
        OptimizerDropDown_2Label        matlab.ui.control.Label
        IterationsforOptimizationEditField  matlab.ui.control.NumericEditField
        IterationsforOptimizationEditFieldLabel  matlab.ui.control.Label
        ConductOptimizationPanel        matlab.ui.container.Panel
        TailoredStartingGuessifavailableCheckBox  matlab.ui.control.CheckBox
        MaterialBudgetEditField         matlab.ui.control.NumericEditField
        MaterialBudgetEditField_2Label  matlab.ui.control.Label
        RunTopologyOptimizationButton   matlab.ui.control.Button
        DebugModeCheckBox               matlab.ui.control.CheckBox
        FilteringProjectionPenaltyPanel  matlab.ui.container.Panel
        MaximumPenaltyofHeavisideProjectionEditField  matlab.ui.control.NumericEditField
        MaximumPenaltyofHeavisideProjectionEditFieldLabel  matlab.ui.control.Label
        RadiusofDensitybasedFilterEditField  matlab.ui.control.NumericEditField
        RadiusofDensitybasedFilterEditField_2Label  matlab.ui.control.Label
        SpecifyDesignDomainPanel        matlab.ui.container.Panel
        VolumeFractionPassiveElementsEditField  matlab.ui.control.NumericEditField
        VolumeFractionPassiveElementsEditFieldLabel  matlab.ui.control.Label
        ClearPassiveElementsButton      matlab.ui.control.Button
        FixedAreaEditField              matlab.ui.control.NumericEditField
        FixedAreaEditFieldLabel         matlab.ui.control.Label
        LoadingAreaEditField            matlab.ui.control.NumericEditField
        LoadingAreaEditField_2Label     matlab.ui.control.Label
        EntireBoundaryEditField         matlab.ui.control.NumericEditField
        EntireBoundaryEditField_2Label  matlab.ui.control.Label
        SetPassiveElementsButton        matlab.ui.control.Button
    end

    
    properties (Access = private)
        MainApp % Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            app.MainApp = mainapp;            
            app.DesignComplianceEditField.Editable = 'off';
            app.DesignVolumeFractionEditField.Editable = 'off';
            app.VolumeFractionPassiveElementsEditField.Editable = 'off';
            app.VolumeFractionPassiveElementsEditField.Value = TopOpti_GetVolumeFractionOfPassiveElements();
            ShowDesignDomain_Public(app.MainApp);  
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            if isvalid(app.MainApp)
                % app.MainApp.SimulationTasksDropDown.Enable = 'on';
                % app.MainApp.FEAwithSolidDesignDomainButton.Enable = 'on';
                MainWindowCtrl(app.MainApp, 1);
                app.MainApp.SimulationTasksDropDown.Value = 'None';
            end
            delete(app)               
        end

        % Button pushed function: SetPassiveElementsButton
        function SetPassiveElementsButtonPushed(app, event)
            global voxelsOnBoundary_;

            app.SpecifyDesignDomainPanel.Enable = 'off';
            app.ConductOptimizationPanel.Enable = 'off';
            pause(1);

            numLayerboundary = app.EntireBoundaryEditField.Value;
            numLayerLoads = app.LoadingAreaEditField.Value;
            numLayerFixation = app.FixedAreaEditField.Value;            
            [voxelsOnBoundary_,~,~] = TopOpti_SetPassiveElements(numLayerboundary, numLayerLoads, numLayerFixation);
            app.VolumeFractionPassiveElementsEditField.Value = TopOpti_GetVolumeFractionOfPassiveElements();
            ShowDesignDomain_Public(app.MainApp);  

            app.SpecifyDesignDomainPanel.Enable = 'on';
            app.ConductOptimizationPanel.Enable = 'on';           
        end

        % Button pushed function: RunTopologyOptimizationButton
        function RunTopologyOptimizationButton_Pushed(app, event)
            global rMin_;
            global pMax_;
            global constraintType_;
            global V_;
            global nLoop_;
            global maxSharpness_;
            global minChange_;
            global optimizer_;
            global move_;
            global complianceDesign_;
            global volumeFractionDesign_;
            global complianceSolid_;
            global DEBUG_;
            global SGopt_;
            global axHandle_;

            app.SpecifyDesignDomainPanel.Enable = 'off';
            app.FilteringProjectionPenaltyPanel.Enable = 'off';
            app.OptimizationProcessSettingsPanel.Enable = 'off';
            app.ConductOptimizationPanel.Enable = 'off';
            app.ResultDisplayPanel.Enable = 'off';
            pause(1);

            GatherLSSandMPsettings(app.MainApp);
            rMin_ = app.RadiusofDensitybasedFilterEditField.Value; rMin_ = max(rMin_, 1.6);
            pMax_ = app.MaximumPenaltyofHeavisideProjectionEditField.Value; pMax_ = max(pMax_, 128);
            constraintType_ = 'Global';
            V_ = app.MaterialBudgetEditField.Value;
            nLoop_ = app.IterationsforOptimizationEditField.Value;
            maxSharpness_ = app.ConvergeSharpnessEditField.Value;
            minChange_ = app.ConvergeMinChangeEditField.Value;
            optimizer_ = app.OptimizerDropDown.Value;
            move_ = app.MovingStepsizeEditField.Value;
            DEBUG_ = app.DebugModeCheckBox.Value;
            SGopt_ = app.TailoredStartingGuessifavailableCheckBox.Value;

            TopOpti_CallTopOpti(axHandle_);            
            ShowDesignDensityLayout_Public(app.MainApp);           
            
            app.SpecifyDesignDomainPanel.Enable = 'on';
            app.FilteringProjectionPenaltyPanel.Enable = 'on';
            app.OptimizationProcessSettingsPanel.Enable = 'on';
            app.ConductOptimizationPanel.Enable = 'on';
            app.ResultDisplayPanel.Enable = 'on';            
                app.DesignComplianceEditField.Value = complianceDesign_;
                app.DesignVolumeFractionEditField.Value = volumeFractionDesign_;
            app.MainApp.ShowComplianceHistoryMenu.Enable = 'on';
            app.MainApp.ShowDesignbyIsosurfaceNotrecommendedMenu.Enable = 'on'; 
            app.MainApp.SolidComplianceEditField.Value = complianceSolid_;            
            app.MainApp.DesignVolEditField.Value = volumeFractionDesign_;
            app.MainApp.DesignComplianceEditField.Value = complianceDesign_;
        end

        % Button pushed function: ClearPassiveElementsButton
        function ClearPassiveElementsButtonPushed(app, event)
            global passiveElements_;
            passiveElements_ = [];
            app.EntireBoundaryEditField.Value = 0;
            app.LoadingAreaEditField.Value = 0;
            app.FixedAreaEditField.Value = 0;
            app.VolumeFractionPassiveElementsEditField.Value = 0.0;
            ShowDesignDomain_Public(app.MainApp);  
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 713 659];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create TopologyOptimizationPanel
            app.TopologyOptimizationPanel = uipanel(app.UIFigure);
            app.TopologyOptimizationPanel.Title = 'Topology Optimization';
            app.TopologyOptimizationPanel.Position = [0 14 700 646];

            % Create SpecifyDesignDomainPanel
            app.SpecifyDesignDomainPanel = uipanel(app.TopologyOptimizationPanel);
            app.SpecifyDesignDomainPanel.Title = 'Specify Design Domain';
            app.SpecifyDesignDomainPanel.Position = [0 515 700 110];

            % Create SetPassiveElementsButton
            app.SetPassiveElementsButton = uibutton(app.SpecifyDesignDomainPanel, 'push');
            app.SetPassiveElementsButton.ButtonPushedFcn = createCallbackFcn(app, @SetPassiveElementsButtonPushed, true);
            app.SetPassiveElementsButton.Position = [548 14 132 23];
            app.SetPassiveElementsButton.Text = 'Set Passive Elements';

            % Create EntireBoundaryEditField_2Label
            app.EntireBoundaryEditField_2Label = uilabel(app.SpecifyDesignDomainPanel);
            app.EntireBoundaryEditField_2Label.HorizontalAlignment = 'right';
            app.EntireBoundaryEditField_2Label.Position = [503 55 91 22];
            app.EntireBoundaryEditField_2Label.Text = 'Entire Boundary';

            % Create EntireBoundaryEditField
            app.EntireBoundaryEditField = uieditfield(app.SpecifyDesignDomainPanel, 'numeric');
            app.EntireBoundaryEditField.ValueDisplayFormat = '%.0f';
            app.EntireBoundaryEditField.Position = [609 55 68 22];

            % Create LoadingAreaEditField_2Label
            app.LoadingAreaEditField_2Label = uilabel(app.SpecifyDesignDomainPanel);
            app.LoadingAreaEditField_2Label.HorizontalAlignment = 'right';
            app.LoadingAreaEditField_2Label.Position = [103 55 76 22];
            app.LoadingAreaEditField_2Label.Text = 'Loading Area';

            % Create LoadingAreaEditField
            app.LoadingAreaEditField = uieditfield(app.SpecifyDesignDomainPanel, 'numeric');
            app.LoadingAreaEditField.ValueDisplayFormat = '%.0f';
            app.LoadingAreaEditField.Position = [194 55 68 22];

            % Create FixedAreaEditFieldLabel
            app.FixedAreaEditFieldLabel = uilabel(app.SpecifyDesignDomainPanel);
            app.FixedAreaEditFieldLabel.HorizontalAlignment = 'right';
            app.FixedAreaEditFieldLabel.Position = [315 55 62 22];
            app.FixedAreaEditFieldLabel.Text = 'Fixed Area';

            % Create FixedAreaEditField
            app.FixedAreaEditField = uieditfield(app.SpecifyDesignDomainPanel, 'numeric');
            app.FixedAreaEditField.ValueDisplayFormat = '%.0f';
            app.FixedAreaEditField.Position = [392 55 68 22];

            % Create ClearPassiveElementsButton
            app.ClearPassiveElementsButton = uibutton(app.SpecifyDesignDomainPanel, 'push');
            app.ClearPassiveElementsButton.ButtonPushedFcn = createCallbackFcn(app, @ClearPassiveElementsButtonPushed, true);
            app.ClearPassiveElementsButton.Position = [389 14 142 23];
            app.ClearPassiveElementsButton.Text = 'Clear Passive Elements';

            % Create VolumeFractionPassiveElementsEditFieldLabel
            app.VolumeFractionPassiveElementsEditFieldLabel = uilabel(app.SpecifyDesignDomainPanel);
            app.VolumeFractionPassiveElementsEditFieldLabel.HorizontalAlignment = 'right';
            app.VolumeFractionPassiveElementsEditFieldLabel.Position = [35 14 199 22];
            app.VolumeFractionPassiveElementsEditFieldLabel.Text = 'Volume Fraction (Passive Elements)';

            % Create VolumeFractionPassiveElementsEditField
            app.VolumeFractionPassiveElementsEditField = uieditfield(app.SpecifyDesignDomainPanel, 'numeric');
            app.VolumeFractionPassiveElementsEditField.Position = [249 14 68 22];

            % Create FilteringProjectionPenaltyPanel
            app.FilteringProjectionPenaltyPanel = uipanel(app.TopologyOptimizationPanel);
            app.FilteringProjectionPenaltyPanel.Title = 'Filtering&Projection&Penalty';
            app.FilteringProjectionPenaltyPanel.Position = [0 413 700 101];

            % Create RadiusofDensitybasedFilterEditField_2Label
            app.RadiusofDensitybasedFilterEditField_2Label = uilabel(app.FilteringProjectionPenaltyPanel);
            app.RadiusofDensitybasedFilterEditField_2Label.HorizontalAlignment = 'right';
            app.RadiusofDensitybasedFilterEditField_2Label.Position = [426 44 166 22];
            app.RadiusofDensitybasedFilterEditField_2Label.Text = 'Radius of Density-based Filter';

            % Create RadiusofDensitybasedFilterEditField
            app.RadiusofDensitybasedFilterEditField = uieditfield(app.FilteringProjectionPenaltyPanel, 'numeric');
            app.RadiusofDensitybasedFilterEditField.Position = [607 44 74 22];
            app.RadiusofDensitybasedFilterEditField.Value = 2.6;

            % Create MaximumPenaltyofHeavisideProjectionEditFieldLabel
            app.MaximumPenaltyofHeavisideProjectionEditFieldLabel = uilabel(app.FilteringProjectionPenaltyPanel);
            app.MaximumPenaltyofHeavisideProjectionEditFieldLabel.HorizontalAlignment = 'right';
            app.MaximumPenaltyofHeavisideProjectionEditFieldLabel.Position = [379 11 228 22];
            app.MaximumPenaltyofHeavisideProjectionEditFieldLabel.Text = 'Maximum Penalty of Heaviside Projection';

            % Create MaximumPenaltyofHeavisideProjectionEditField
            app.MaximumPenaltyofHeavisideProjectionEditField = uieditfield(app.FilteringProjectionPenaltyPanel, 'numeric');
            app.MaximumPenaltyofHeavisideProjectionEditField.ValueDisplayFormat = '%.0f';
            app.MaximumPenaltyofHeavisideProjectionEditField.Position = [622 11 60 22];
            app.MaximumPenaltyofHeavisideProjectionEditField.Value = 128;

            % Create ConductOptimizationPanel
            app.ConductOptimizationPanel = uipanel(app.TopologyOptimizationPanel);
            app.ConductOptimizationPanel.Title = 'Conduct Optimization';
            app.ConductOptimizationPanel.Position = [0 82 700 111];

            % Create DebugModeCheckBox
            app.DebugModeCheckBox = uicheckbox(app.ConductOptimizationPanel);
            app.DebugModeCheckBox.Text = 'Debug Mode';
            app.DebugModeCheckBox.Position = [381 17 91 22];

            % Create RunTopologyOptimizationButton
            app.RunTopologyOptimizationButton = uibutton(app.ConductOptimizationPanel, 'push');
            app.RunTopologyOptimizationButton.ButtonPushedFcn = createCallbackFcn(app, @RunTopologyOptimizationButton_Pushed, true);
            app.RunTopologyOptimizationButton.FontWeight = 'bold';
            app.RunTopologyOptimizationButton.Position = [516 17 171 23];
            app.RunTopologyOptimizationButton.Text = 'Run Topology Optimization';

            % Create MaterialBudgetEditField_2Label
            app.MaterialBudgetEditField_2Label = uilabel(app.ConductOptimizationPanel);
            app.MaterialBudgetEditField_2Label.HorizontalAlignment = 'right';
            app.MaterialBudgetEditField_2Label.Position = [481 52 89 22];
            app.MaterialBudgetEditField_2Label.Text = 'Material Budget';

            % Create MaterialBudgetEditField
            app.MaterialBudgetEditField = uieditfield(app.ConductOptimizationPanel, 'numeric');
            app.MaterialBudgetEditField.Position = [585 52 100 22];
            app.MaterialBudgetEditField.Value = 0.4;

            % Create TailoredStartingGuessifavailableCheckBox
            app.TailoredStartingGuessifavailableCheckBox = uicheckbox(app.ConductOptimizationPanel);
            app.TailoredStartingGuessifavailableCheckBox.Text = 'Tailored Starting Guess (if available)';
            app.TailoredStartingGuessifavailableCheckBox.Position = [103 17 216 22];

            % Create OptimizationProcessSettingsPanel
            app.OptimizationProcessSettingsPanel = uipanel(app.TopologyOptimizationPanel);
            app.OptimizationProcessSettingsPanel.Title = 'Optimization Process Settings';
            app.OptimizationProcessSettingsPanel.Position = [0 192 700 222];

            % Create IterationsforOptimizationEditFieldLabel
            app.IterationsforOptimizationEditFieldLabel = uilabel(app.OptimizationProcessSettingsPanel);
            app.IterationsforOptimizationEditFieldLabel.HorizontalAlignment = 'right';
            app.IterationsforOptimizationEditFieldLabel.Position = [420 171 148 22];
            app.IterationsforOptimizationEditFieldLabel.Text = '#Iterations for Optimization';

            % Create IterationsforOptimizationEditField
            app.IterationsforOptimizationEditField = uieditfield(app.OptimizationProcessSettingsPanel, 'numeric');
            app.IterationsforOptimizationEditField.ValueDisplayFormat = '%.0f';
            app.IterationsforOptimizationEditField.Position = [583 171 100 22];
            app.IterationsforOptimizationEditField.Value = 50;

            % Create OptimizerDropDown_2Label
            app.OptimizerDropDown_2Label = uilabel(app.OptimizationProcessSettingsPanel);
            app.OptimizerDropDown_2Label.HorizontalAlignment = 'right';
            app.OptimizerDropDown_2Label.Position = [513 58 56 22];
            app.OptimizerDropDown_2Label.Text = 'Optimizer';

            % Create OptimizerDropDown
            app.OptimizerDropDown = uidropdown(app.OptimizationProcessSettingsPanel);
            app.OptimizerDropDown.Items = {'OC', 'MMA'};
            app.OptimizerDropDown.Position = [584 58 100 22];
            app.OptimizerDropDown.Value = 'OC';

            % Create MovingStepsizeEditField_2Label
            app.MovingStepsizeEditField_2Label = uilabel(app.OptimizationProcessSettingsPanel);
            app.MovingStepsizeEditField_2Label.HorizontalAlignment = 'right';
            app.MovingStepsizeEditField_2Label.Position = [476 17 93 22];
            app.MovingStepsizeEditField_2Label.Text = 'Moving Stepsize';

            % Create MovingStepsizeEditField
            app.MovingStepsizeEditField = uieditfield(app.OptimizationProcessSettingsPanel, 'numeric');
            app.MovingStepsizeEditField.Position = [584 17 100 22];
            app.MovingStepsizeEditField.Value = 0.1;

            % Create ConvergeSharpnessEditFieldLabel
            app.ConvergeSharpnessEditFieldLabel = uilabel(app.OptimizationProcessSettingsPanel);
            app.ConvergeSharpnessEditFieldLabel.HorizontalAlignment = 'right';
            app.ConvergeSharpnessEditFieldLabel.Position = [499 133 118 22];
            app.ConvergeSharpnessEditFieldLabel.Text = 'Converge Sharpness';

            % Create ConvergeSharpnessEditField
            app.ConvergeSharpnessEditField = uieditfield(app.OptimizationProcessSettingsPanel, 'numeric');
            app.ConvergeSharpnessEditField.Position = [632 133 51 22];
            app.ConvergeSharpnessEditField.Value = 0.05;

            % Create ConvergeMinChangeEditFieldLabel
            app.ConvergeMinChangeEditFieldLabel = uilabel(app.OptimizationProcessSettingsPanel);
            app.ConvergeMinChangeEditFieldLabel.HorizontalAlignment = 'right';
            app.ConvergeMinChangeEditFieldLabel.Position = [496 95 125 22];
            app.ConvergeMinChangeEditFieldLabel.Text = 'Converge Min Change';

            % Create ConvergeMinChangeEditField
            app.ConvergeMinChangeEditField = uieditfield(app.OptimizationProcessSettingsPanel, 'numeric');
            app.ConvergeMinChangeEditField.Position = [636 95 47 22];
            app.ConvergeMinChangeEditField.Value = 0.0001;

            % Create ResultDisplayPanel
            app.ResultDisplayPanel = uipanel(app.TopologyOptimizationPanel);
            app.ResultDisplayPanel.Title = 'Result Display';
            app.ResultDisplayPanel.Position = [0 1 700 82];

            % Create DesignComplianceEditFieldLabel
            app.DesignComplianceEditFieldLabel = uilabel(app.ResultDisplayPanel);
            app.DesignComplianceEditFieldLabel.HorizontalAlignment = 'right';
            app.DesignComplianceEditFieldLabel.Position = [79 16 109 22];
            app.DesignComplianceEditFieldLabel.Text = 'Design Compliance';

            % Create DesignComplianceEditField
            app.DesignComplianceEditField = uieditfield(app.ResultDisplayPanel, 'numeric');
            app.DesignComplianceEditField.Position = [202 16 100 22];

            % Create DesignVolumeFractionEditFieldLabel
            app.DesignVolumeFractionEditFieldLabel = uilabel(app.ResultDisplayPanel);
            app.DesignVolumeFractionEditFieldLabel.HorizontalAlignment = 'right';
            app.DesignVolumeFractionEditFieldLabel.Position = [440 16 132 22];
            app.DesignVolumeFractionEditFieldLabel.Text = 'Design Volume Fraction';

            % Create DesignVolumeFractionEditField
            app.DesignVolumeFractionEditField = uieditfield(app.ResultDisplayPanel, 'numeric');
            app.DesignVolumeFractionEditField.Position = [586 16 100 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Mtd_TopologyOptimization(varargin)

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