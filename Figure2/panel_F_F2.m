function panel_F_F2(decoding,confusion,accuracy,xrange,ax)

    stimlabel = {0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1};
    n_stim = size(decoding{1}{1},1);
    n_trials = size(decoding{1}{1},2);
    cspan = [0 50];
    col = [0 0.447 0.741; 0.85 0.325 0.098];

    for cond = 1:2
        axes(ax(cond));
        C = cat(3,confusion{cond}{end,:,:});
        C_avg = 100*nanmean(C,3)/n_trials;
        imagesc(C_avg);
        axis square;
        hold on;
        for i = 1:length(xrange)
            for j = 1:length(xrange)
                line([0.5, i+0.5], [j+0.5, j+0.5],...
                    'Color', 'k',...
                    'LineWidth', 1);
                line([i+0.5, i+0.5], [0.5, j+0.5],...
                    'Color', 'k',...
                    'LineWidth', 1);
            end
        end
        
        if cond==2
            c = mycolorbar('%');
            c.Ticks = [0 25 50];
            % manually place the colorbar because Matlab is broken with axis square
            c.Position = [0.89 0.26 0.025 0.705];
            c.Label.Position(1) = 0;
            c.Label.Position(2) = 25;
            c.Label.Color = 'w';
            c.Label.FontWeight = 'bold';
            c.Label.FontSize = 8;
            ax(cond).YTick = 1:n_stim;
            ax(cond).YTickLabel = {};
        else
            ax(cond).YTick = 1:n_stim;
            ax(cond).YTickLabel = stimlabel(xrange);
            ylabel('Stim. presented');
        end
        colormap(1-gray);
        ax(cond).XTick = 1:n_stim;
        ax(cond).XTickLabel =  stimlabel(xrange);
        ax(cond).XTickLabelRotation = 45;
        xlabel('Stim. classified as');
        caxis(cspan);
        if cond == 1
            title(ax(cond),'wild type','Color',col(1,:));
        else
            title(ax(cond),'{\it{}mecp2}-mutant','Color',col(2,:));
        end
    end
end