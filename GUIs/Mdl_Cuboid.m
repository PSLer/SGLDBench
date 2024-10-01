classdef Mdl_Cuboid < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        CreateCuboidDesignDomainPanel  matlab.ui.container.Panel
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
            if 0>=xDim || 0>=yDim || 0 >= zDim, return; end
            Shape_BuiltInCuboid(xDim, yDim, zDim);
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
            app.UIFigure.Position = [100 100 414 252];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create CreateCuboidDesignDomainPanel
            app.CreateCuboidDesignDomainPanel = uipanel(app.UIFigure);
            app.CreateCuboidDesignDomainPanel.Title = 'Create Cuboid Design Domain';
            app.CreateCuboidDesignDomainPanel.Position = [2 -2 416 255];

            % Create DimensionXEditFieldLabel
            app.DimensionXEditFieldLabel = uilabel(app.CreateCuboidDesignDomainPanel);
            app.DimensionXEditFieldLabel.HorizontalAlignment = 'right';
            app.DimensionXEditFieldLabel.Position = [211 179 73 22];
            app.DimensionXEditFieldLabel.Text = 'Dimension X';

            % Create DimensionXEditField
            app.DimensionXEditField = uieditfield(app.CreateCuboidDesignDomainPanel, 'numeric');
            app.DimensionXEditField.Position = [299 179 100 22];
            app.DimensionXEditField.Value = 2;

            % Create DimensionYEditFieldLabel
            app.DimensionYEditFieldLabel = uilabel(app.CreateCuboidDesignDomainPanel);
            app.DimensionYEditFieldLabel.HorizontalAlignment = 'right';
            app.DimensionYEditFieldLabel.Position = [212 138 73 22];
            app.DimensionYEditFieldLabel.Text = 'Dimension Y';

            % Create DimensionYEditField
            app.DimensionYEditField = uieditfield(app.CreateCuboidDesignDomainPanel, 'numeric');
            app.DimensionYEditField.Position = [300 138 100 22];
            app.DimensionYEditField.Value = 1;

            % Create DimensionZEditFieldLabel
            app.DimensionZEditFieldLabel = uilabel(app.CreateCuboidDesignDomainPanel);
            app.DimensionZEditFieldLabel.HorizontalAlignment = 'right';
            app.DimensionZEditFieldLabel.Position = [212 97 72 22];
            app.DimensionZEditFieldLabel.Text = 'Dimension Z';

            % Create DimensionZEditField
            app.DimensionZEditField = uieditfield(app.CreateCuboidDesignDomainPanel, 'numeric');
            app.DimensionZEditField.Position = [299 97 100 22];
            app.DimensionZEditField.Value = 1;

            % Create CreateButton
            app.CreateButton = uibutton(app.CreateCuboidDesignDomainPanel, 'push');
            app.CreateButton.ButtonPushedFcn = createCallbackFcn(app, @CreateButtonPushed, true);
            app.CreateButton.Position = [301 51 100 23];
            app.CreateButton.Text = 'Create';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Mdl_Cuboid(varargin)

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