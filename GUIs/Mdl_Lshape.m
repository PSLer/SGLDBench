classdef Mdl_Lshape < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        CreateLshapeDesignDomainPanel  matlab.ui.container.Panel
        ControlPointZEditField         matlab.ui.control.NumericEditField
        ControlPointZEditFieldLabel    matlab.ui.control.Label
        ControlPointXEditField         matlab.ui.control.NumericEditField
        ControlPointXEditFieldLabel    matlab.ui.control.Label
        CreateButton                   matlab.ui.control.Button
        DimensionZEditField            matlab.ui.control.NumericEditField
        DimensionZEditFieldLabel       matlab.ui.control.Label
        DimensionYEditField            matlab.ui.control.NumericEditField
        DimensionYEditFieldLabel       matlab.ui.control.Label
        DimensionXEditField            matlab.ui.control.NumericEditField
        DimensionXEditFieldLabel       matlab.ui.control.Label
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
            xDim = app.DimensionXEditField.Value;
            yDim = app.DimensionYEditField.Value;
            zDim = app.DimensionZEditField.Value;
            cPointX = app.ControlPointXEditField.Value;
            cPointZ = app.ControlPointZEditField.Value;
            if 0>=xDim || 0>=yDim || 0 >= zDim || 0 >= cPointX || 0 >= cPointZ, return; end
            if cPointX>=xDim || cPointZ >=zDim, return; end
            Shape_BuiltInLshape(xDim, yDim, zDim, cPointX, cPointZ);
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
            app.UIFigure.Position = [100 100 478 279];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create CreateLshapeDesignDomainPanel
            app.CreateLshapeDesignDomainPanel = uipanel(app.UIFigure);
            app.CreateLshapeDesignDomainPanel.Title = 'Create L-shape Design Domain';
            app.CreateLshapeDesignDomainPanel.Position = [2 0 478 280];

            % Create DimensionXEditFieldLabel
            app.DimensionXEditFieldLabel = uilabel(app.CreateLshapeDesignDomainPanel);
            app.DimensionXEditFieldLabel.HorizontalAlignment = 'right';
            app.DimensionXEditFieldLabel.Position = [266 204 73 22];
            app.DimensionXEditFieldLabel.Text = 'Dimension X';

            % Create DimensionXEditField
            app.DimensionXEditField = uieditfield(app.CreateLshapeDesignDomainPanel, 'numeric');
            app.DimensionXEditField.Position = [354 204 100 22];
            app.DimensionXEditField.Value = 100;

            % Create DimensionYEditFieldLabel
            app.DimensionYEditFieldLabel = uilabel(app.CreateLshapeDesignDomainPanel);
            app.DimensionYEditFieldLabel.HorizontalAlignment = 'right';
            app.DimensionYEditFieldLabel.Position = [267 163 73 22];
            app.DimensionYEditFieldLabel.Text = 'Dimension Y';

            % Create DimensionYEditField
            app.DimensionYEditField = uieditfield(app.CreateLshapeDesignDomainPanel, 'numeric');
            app.DimensionYEditField.Position = [355 163 100 22];
            app.DimensionYEditField.Value = 50;

            % Create DimensionZEditFieldLabel
            app.DimensionZEditFieldLabel = uilabel(app.CreateLshapeDesignDomainPanel);
            app.DimensionZEditFieldLabel.HorizontalAlignment = 'right';
            app.DimensionZEditFieldLabel.Position = [267 122 72 22];
            app.DimensionZEditFieldLabel.Text = 'Dimension Z';

            % Create DimensionZEditField
            app.DimensionZEditField = uieditfield(app.CreateLshapeDesignDomainPanel, 'numeric');
            app.DimensionZEditField.Position = [354 122 100 22];
            app.DimensionZEditField.Value = 100;

            % Create CreateButton
            app.CreateButton = uibutton(app.CreateLshapeDesignDomainPanel, 'push');
            app.CreateButton.ButtonPushedFcn = createCallbackFcn(app, @CreateButtonPushed, true);
            app.CreateButton.Position = [357 29 100 23];
            app.CreateButton.Text = 'Create';

            % Create ControlPointXEditFieldLabel
            app.ControlPointXEditFieldLabel = uilabel(app.CreateLshapeDesignDomainPanel);
            app.ControlPointXEditFieldLabel.HorizontalAlignment = 'right';
            app.ControlPointXEditFieldLabel.Position = [33 76 86 22];
            app.ControlPointXEditFieldLabel.Text = 'Control Point X';

            % Create ControlPointXEditField
            app.ControlPointXEditField = uieditfield(app.CreateLshapeDesignDomainPanel, 'numeric');
            app.ControlPointXEditField.Position = [134 76 100 22];
            app.ControlPointXEditField.Value = 50;

            % Create ControlPointZEditFieldLabel
            app.ControlPointZEditFieldLabel = uilabel(app.CreateLshapeDesignDomainPanel);
            app.ControlPointZEditFieldLabel.HorizontalAlignment = 'right';
            app.ControlPointZEditFieldLabel.Position = [256 76 85 22];
            app.ControlPointZEditFieldLabel.Text = 'Control Point Z';

            % Create ControlPointZEditField
            app.ControlPointZEditField = uieditfield(app.CreateLshapeDesignDomainPanel, 'numeric');
            app.ControlPointZEditField.Position = [356 76 100 22];
            app.ControlPointZEditField.Value = 50;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Mdl_Lshape(varargin)

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