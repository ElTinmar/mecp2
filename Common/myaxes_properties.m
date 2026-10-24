function myaxes_properties(ax_list)
    for ax = ax_list
        ax.LabelFontSizeMultiplier = 1;
        ax.TitleFontSizeMultiplier = 1;
        ax.Box = 'on';
        ax.FontSize = 8;
        ax.XLabel.FontSize = 10;
        ax.YLabel.FontSize = 10;
        ax.Title.FontSize = 10;
        ax.LineWidth = 1;
        ax.FontName = 'Arial';
        
        %% undocumented properties
        ax.XRuler.TickLabelGapOffset = 0.1;
        ax.YRuler.TickLabelGapOffset = 0.1;
        ax.LooseInset = max(ax.TightInset, 0.02);

        % make colorbar same size as plot
        if ~isempty(ax.Colorbar)
            ax.Colorbar.FontSize = 8;
        end        
    end
end

