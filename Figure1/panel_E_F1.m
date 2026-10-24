function panel_E_F1(activity,corrdata,ax)

    yspan = [1e-7 1e-1];
    xspan= [-.1 .85];
    xticks = [0 .2 .4 .6 .8];
    col = [0 0.447 0.741; 0.85 0.325 0.098];

    plot_corr_distrib(1)
    plot_corr_distrib(2)

    function plot_corr_distrib(i)
        axes(ax(i))
        hold on

        if i==1
            distrib = mean(corrdata.wt.corr_dist);
            distrib_null = mean(reshape(corrdata.wt.corr_dist_null,[],1000));
        else
            distrib = mean(corrdata.mecp2.corr_dist);
            distrib_null = mean(reshape(corrdata.mecp2.corr_dist_null,[],1000));
        end

        % plot data
        p = semilogy(corrdata.corr_centers,distrib);
        p.Marker = '.';
        p.Color = col(i,:);
        p.LineStyle = 'none';
        p.MarkerSize = 5;

        % plot null model
        p0 = semilogy(corrdata.corr_centers,distrib_null);
        p0.Marker = '.';
        p0.Color = [.2 .2 .2];
        p0.LineStyle = 'none';
        p0.MarkerSize = 5;

        p1 = semilogy([0 0],yspan,'k');

        ylim(yspan);
        xlim(xspan);
        xlabel('Correlation');

        if i==1
            ylabel('<Probability>_{fish}');
            ax(i).YTick = [1e-6 1e-4 1e-2 1e0];
        else
            ax(i).YTickLabels = {};
        end
        ax(i).YScale = 'log';
        ax(i).XTick = xticks;

        if i==1
            %lgd = legend([p p0],{'wild type','surrogate'},'FontSize',8,'Box','off');
            %lgd.Position(1) = 0.225;
            text(0.8, 0.02, 'wild type', 'FontSize', 8, 'FontWeight', 'bold', 'Color', col(i,:),'HorizontalAlignment', 'right');
            text(0.8, 0.005, 'surrogate', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'right');
        else
            %lgd = legend([p p0],{'{\it{}mecp2}-mutant','surrogate'},'FontSize',8,'Box','off');
            %lgd.Position(1) = 0.625;
            text(0.8, 0.02, '{\it{}mecp2}-mutant', 'FontSize', 8, 'FontWeight', 'bold', 'Color', col(i,:), 'HorizontalAlignment', 'right');
            text(0.8, 0.005, 'surrogate', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'right');
        end
        %lgd.ItemTokenSize = [5 10];
        
    end
end