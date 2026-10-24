function panel_E_F3(decoding,confusion,accuracy,stimlabel,xrange,ax)

    n_stim = size(decoding{1}{1},1);
    n_trials = size(decoding{1}{1},2);
    cspan = [0 60];
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
            c.Ticks = [0 30 60];
            % manually place the colorbar because Matlab is broken with axis square
            c.Position = [0.92 0.262 0.025 0.701];
            c.Label.Position(1) = 0;
            c.Label.Position(2) = 30;
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
        xlabel('Stim. classified as');
        caxis(cspan);
        if cond == 1
            title(ax(cond),'wild type','Color',col(1,:));
        else
            title(ax(cond),'{\it{}mecp2}-mutant','Color',col(2,:));
        end
    end
end