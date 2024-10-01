classdef Mdl_Cylinder < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        CreateCylinderShellDesignDomainPanel  matlab.ui.container.Panel
        HeightEditField                 matlab.ui.control.NumericEditField
        HeightEditFieldLabel            matlab.ui.control.Label
        OuterRadiusEditField            matlab.ui.control.NumericEditField
        OuterRadiusEditFieldLabel       matlab.ui.control.Label
        CreateButton                    matlab.ui.control.Button
        InnerRadiusRatioEditField       matlab.ui.control.NumericEditField
        InnerRadiusRatioEditFieldLabel  matlab.ui.control.Label
    end

    
    properties (Access = private)
        MainApp % Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            app.MainApp = mainapp;
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            if isvalid(app.MainApp)
                app.MainApp.ImportMenu.Enable = 'on';
                app.MainApp.VisualizationMenu.Enable = 'on';
                app.MainApp.DomainVoxelizationPanel.Enable = 'on';                 
            end
            delete(app)            
        end

        % Button pushed function: CreateButton
        function CreateButtonPushed(app, event)
            outerRadius = app.OuterRadiusEditField.Value;
            height = app.HeightEditField.Value;
            ratioInnerRadius = app.InnerRadiusRatioEditField.Value;
            if ratioInnerRadius<0 || ratioInnerRadius>1, return; end
            if outerRadius<=0, return; end
            if height<=0, return; end
            Shape_BuiltInCylinder(outerRadius, ratioInnerRadius, height);
            app.MainApp.ShowTriSurfaceMesh_public();
            % app.MainApp.BuiltinShapesMenu.Enable = 'on';
            % delete(app);
            app.MainApp.VoxelizingButton.Enable = 'on';
            app.MainApp.ShowInputTriangularSurfaceMeshMenu.Enable = 'on';            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 433 225];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create CreateCylinderShellDesignDomainPanel
            app.CreateCylinderShellDesignDomainPanel = uipanel(app.UIFigure);
            app.CreateCylinderShellDesignDomainPanel.Title = 'Create Cylinder (Shell) Design Domain';
            app.CreateCylinderShellDesignDomainPanel.Position = [2 -1 434 227];

            % Create InnerRadiusRatioEditFieldLabel
            app.InnerRadiusRatioEditFieldLabel = uilabel(app.CreateCylinderShellDesignDomainPanel);
            app.InnerRadiusRatioEditFieldLabel.HorizontalAlignment = 'right';
            app.InnerRadiusRatioEditFieldLabel.Position = [196 77 104 22];
            app.InnerRadiusRatioEditFieldLabel.Text = 'Inner Radius Ratio';

            % Create InnerRadiusRatioEditField
            app.InnerRadiusRatioEditField = uieditfield(app.CreateCylinderShellDesignDomainPanel, 'numeric');
            app.InnerRadiusRatioEditField.Position = [315 77 100 22];
            app.InnerRadiusRatioEditField.Value = 0.8;

            % Create CreateButton
            app.CreateButton = uibutton(app.CreateCylinderShellDesignDomainPanel, 'push');
            app.CreateButton.ButtonPushedFcn = createCallbackFcn(app, @CreateButtonPushed, true);
            app.CreateButton.Position = [314 34 100 23];
            app.CreateButton.Text = 'Create';

            % Create OuterRadiusEditFieldLabel
            app.OuterRadiusEditFieldLabel = uilabel(app.CreateCylinderShellDesignDomainPanel);
            app.OuterRadiusEditFieldLabel.HorizontalAlignment = 'right';
            app.OuterRadiusEditFieldLabel.Position = [223 157 76 22];
            app.OuterRadiusEditFieldLabel.Text = 'Outer Radius';

            % Create OuterRadiusEditField
            app.OuterRadiusEditField = uieditfield(app.CreateCylinderShellDesignDomainPanel, 'numeric');
            app.OuterRadiusEditField.Position = [314 157 100 22];
            app.OuterRadiusEditField.Value = 1;

            % Create HeightEditFieldLabel
            app.HeightEditFieldLabel = uilabel(app.CreateCylinderShellDesignDomainPanel);
            app.HeightEditFieldLabel.HorizontalAlignment = 'right';
            app.HeightEditFieldLabel.Position = [260 116 40 22];
            app.HeightEditFieldLabel.Text = 'Height';

            % Create HeightEditField
            app.HeightEditField = uieditfield(app.CreateCylinderShellDesignDomainPanel, 'numeric');
            app.HeightEditField.Position = [315 116 100 22];
            app.HeightEditField.Value = 5;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Mdl_Cylinder(varargin)

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