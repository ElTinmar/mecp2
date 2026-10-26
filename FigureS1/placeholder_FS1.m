function placeholder_FS1(ax,msg)
    axes(ax)
    text(0.5,0.5,msg,'HorizontalAlignment','center','Color',[0.5 0.5 0.5],'FontSize',8)
    ax.XTick = [];
    ax.YTick = [];
    ax.Box = 'on';
    ax.XColor = [0.7 0.7 0.7];
    ax.YColor = [0.7 0.7 0.7];
    xlim([0 1]); ylim([0 1]);
end
