function f = myfigure()

    f = figure();
    pause(0.1)
    f.PaperUnits = 'centimeters';
    f.PaperSize = [21 29.7];
    f.PaperPosition = [0 0 21 29.7];
    f.PaperPositionMode = 'manual';
    f.PaperType = 'A4';
    f.Units = 'centimeters';
    f.Resize = 'off';
    f.Position = [0 0 21 29.7];
    f.MenuBar = 'none';
    f.ToolBar = 'none';
    f.Color = 'w';
    f.InvertHardcopy = 'off';
end
