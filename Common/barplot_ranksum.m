function [p,mx,sx,my,sy] = barplot_ranksum(ax,x,y,ylab,xticklab,col,varargin)

% make sure x and y are row vectors ---------------------------------------
if iscolumn(x), x = x'; end
if iscolumn(y), y = y'; end

% ranksum -----------------------------------------------------------------
[p,h,stats] = ranksum(x,y,varargin{:});

% means and standard dev --------------------------------------------------
mx = nanmean(x);
my = nanmean(y);
sx = nanstd(x)/sqrt(sum(~isnan(x)));
sy = nanstd(y)/sqrt(sum(~isnan(y)));
swarm_x = [ones(1,length(x)) 2*ones(1,length(y))];
swarm_y = [x, y];

% constants ---------------------------------------------------------------
linecol = [0 0 0];
lwidth = 1;
prct_offset = 0.05;

% plot figure -------------------------------------------------------------
axes(ax)
ax.XTick = [1 2];
hold on
my_bar(1,mx,col(1,:),linecol);
my_bar(2,my,col(2,:),linecol);
my_errorbar(1,mx,sx,linecol);
my_errorbar(2,my,sy,linecol);
comparison(mx,sx,my,sy,p,linecol);
swarmchart(swarm_x,swarm_y,'.k','XJitterWidth',0.5);
lims = ylim();
new_lim = max(lims(2), max(swarm_y));
ylim([lims(1) new_lim])
xlim([0.5 2.5])
ylabel(ylab)
ax.XTickLabel = xticklab;

    function my_bar(x,m,facecol,edgecol)
        bar(x,m,...
            'FaceColor',facecol,...
            'EdgeColor',edgecol,...
            'LineWidth',lwidth);
    end

    function  my_errorbar(x,m,s,edgecol)
        plot([x x],[m m+s],...
            'Color',edgecol,...
            'LineWidth',lwidth);
        plot([x-0.05 x+0.05],[m+s m+s],...
            'Color',edgecol,...
            'LineWidth',lwidth);
    end

    function comparison(mx,sx,my,sy,p,edgecol)

        lims = ylim();
        height = diff(lims);
        offset = prct_offset * height;
        top = max(mx+sx ,my+sy) + 2*offset;
        bot1 = mx+sx+offset;
        bot2 = my+sy+offset;
        xp = 1.5;
        yp = top+offset;

        % plot bridge
        plot([1 1],[bot1 top],...
            'Color',edgecol,...
            'LineWidth',lwidth);
        plot([1 2],[top top],...
            'Color',edgecol,...
            'LineWidth',lwidth);
        plot([2 2],[top bot2],...
            'Color',edgecol,...
            'LineWidth',lwidth);
        
        % asterisk p-val
        if p < 0.001
            my_asterisk(xp,yp,'***',12)
        elseif p < 0.01
            my_asterisk(xp,yp,'**',12)
        elseif p < 0.05
            my_asterisk(xp,yp,'*',12)
        else
            my_asterisk(xp,yp+offset,'ns',8)
        end
    
        ylim([lims(1) yp+3*offset])

        function my_asterisk(x,y,txt,font_size)
            text(x,y,txt,...
                'HorizontalAlignment','center',...
                'Color',edgecol,...
                'FontSize',font_size)
        end
    end
end