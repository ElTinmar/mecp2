function c = mycolorbar(label)
    c = colorbar('eastoutside');
    ylabel(c,label);
    c.FontSize = 10;
    c.LineWidth = 1;

    c.Units = 'centimeters';
    c.Position(3) = 0.3;
    c.Units = 'normalized';
end