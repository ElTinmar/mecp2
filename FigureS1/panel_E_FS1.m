function panel_E_FS1(ax)
%%% rostro-caudal distribution: needs the per-cell distance to the tectal midline (VerTodo files)
    TYPES = ["vglut2a","gad1b"];
    for k = 1:2
        axes(ax(k))
        text(0.5,0.5,{'needs distance to midline','(not in ALL.mat)'},...
            'HorizontalAlignment','center','Color',[0.5 0.5 0.5],'FontSize',8)
        xlim([0 1])
        ylabel(sprintf('# {\\it{}%s} cells',TYPES(k)))
        if k == 2, xlabel('Rostro-caudal axis (norm.)'), else, ax(k).XTickLabel = []; end
        ax(k).YTick = [];
        ax(k).PlotBoxAspectRatio = [1.3 1 1];
    end
end
