function [p_local,p_global,anova_table,p_raw] = property_perStim_ranksum(...
    ax,...
    quantity,...
    yspan,...
    ylab,...
    xticklab,...
    xrange,...
    xlab,...
    col)

    quantity_avg{1} = nanmean(quantity{1});
    quantity_avg{2} = nanmean(quantity{2});
    quantity_std{1} = nanstd(quantity{1})./sqrt(sum(~isnan(quantity{1})));
    quantity_std{2} = nanstd(quantity{2})./sqrt(sum(~isnan(quantity{2})));
    
    p_raw = stats_local();
    p_local = NaN(size(p_raw));
    p_local(xrange) = holm(p_raw(xrange));
    [p_global, anova_table] = mixed_anova(quantity{1}(:,xrange),quantity{2}(:,xrange));
    plot_results()

    function plot_results()
        axes(ax)
        for i = 1:2
            y1 = quantity_avg{i}-quantity_std{i};
            y2 = quantity_avg{i}+quantity_std{i};
            plot(xrange,quantity_avg{i}(xrange),'color',col(i,:),'LineWidth',2);
            hold on
            X = [xrange,fliplr(xrange)];
            Y = [y1(xrange),fliplr(y2(xrange))];
            torm = isnan(Y);
            X(torm) = [];
            Y(torm) = [];
            h1 = fill(X,Y,col(i,:),'EdgeColor','none');
            set(h1,'FaceAlpha',0.5);
        end
        ax = gca();
        ax.XTick = xrange(1:2:end);
        ax.XTickLabelRotation = 0;
        ax.XTickLabel = xticklab(xrange(1:2:end));
        xlim([xrange(1)-0.6 xrange(end)+0.6])
        ylim(yspan)
        ylabel(ylab)
        xlabel(xlab)
        axis square
        
        % add p-values to the plot
        for i = find(p_local <= 0.05)
            X = i + [-0.5 0.5 0.5 -0.5];
            Y = [yspan(1) yspan(1) yspan(2) yspan(2)];
            h = fill(X,Y,[0.3 0.3 0.3],'EdgeColor','none');
            h.FaceAlpha = 0.3;
        end
    end

    function p = stats_local()
        % local ranksum test
        p = [];
        for stim = 1:size(quantity{1},2)
            y1 = quantity{1}(:,stim);
            y2 = quantity{2}(:,stim);
            if all(isnan(y1)) | all(isnan(y2))
                p(stim) = NaN;
            else
                p(stim) = ranksum(y1,y2);
            end
        end
    end

end