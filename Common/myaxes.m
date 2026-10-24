function ax = myaxes(parent,pos)
    ax = axes(parent);
    ax.Units = 'normalized';
    ax.PositionConstraint = 'innerposition';
    ax.InnerPosition = pos;
    ax.PlotBoxAspectRatio = [1 1 1];
end