function distribution_FS1(cells,var,range,n_bins,scale,xlab,ax)
%%% distribution_FS1(cells,var,range,n_bins,scale,xlab,ax)
% number of vglut2a (ax(1)) and gad1b (ax(2)) cells per bin of cells.(var), mean +- SEM across larvae.
% Statistics on N_STAT coarser bins: mixed-design ANOVA genotype x bin (the interaction tests the
% shape of the distribution) and Holm-corrected rank-sum tests per bin (shaded).

    N_STAT = 10;
    col = [0 0.447 0.741; 0.85 0.325 0.098];
    TYPES = ["vglut2a","gad1b"];
    larvae = unique(cells.larva);
    wt = arrayfun(@(l) cells.genotype(find(cells.larva == l,1)) == "wt",larvae);
    edges = linspace(range(1),range(2),n_bins+1);
    stat_edges = linspace(range(1),range(2),N_STAT+1);
    x = scale * (edges(1:end-1) + edges(2:end))/2;

    for k = 1:2
        H = zeros(numel(larvae),n_bins);
        Hs = zeros(numel(larvae),N_STAT);
        for l = 1:numel(larvae)
            v = cells.(var)(cells.larva == larvae(l) & cells.type == TYPES(k));
            H(l,:) = histcounts(v,edges);
            Hs(l,:) = histcounts(v,stat_edges);
        end

        [p_anova,tbl] = mixed_anova(Hs(wt,:),Hs(~wt,:));
        p_bin = arrayfun(@(b) ranksum(Hs(wt,b),Hs(~wt,b)),1:N_STAT);
        p_holm = holm(p_bin);
        fprintf('%s, cells per bin (%d bins), mixed-design ANOVA (within effects GG-corrected):\n',TYPES(k),N_STAT);
        fprintf('    MeCP2/WT: %.3g\n    Bin: %.3g\n    Interaction (distribution shape): %.3g\n',p_anova);
        disp(tbl.Properties.Description)
        disp(tbl)
        fprintf('Post-hoc: bilateral ranksum test per bin (edges %s):\n',mat2str(scale*stat_edges,3));
        disp(p_bin)
        disp('Holm-corrected:')
        disp(p_holm)
        disp(' ')

        axes(ax(k))
        hold on
        yhi = 1.15 * max(mean(H) + std(H)/sqrt(size(H,1)),[],'all');
        for b = find(p_holm < 0.05)
            fill(scale*stat_edges([b b+1 b+1 b]),[0 0 yhi yhi],[0.3 0.3 0.3],'EdgeColor','none','FaceAlpha',0.3);
        end
        for g = 1:2
            Y = H(wt == (g == 1),:);
            m = mean(Y);
            se = std(Y)/sqrt(size(Y,1));
            fill([x fliplr(x)],[m-se fliplr(m+se)],col(g,:),'EdgeColor','none','FaceAlpha',0.3);
            plot(x,m,'Color',col(g,:),'LineWidth',1.5);
        end
        xlim(scale*range)
        ylim([0 yhi])
        ylabel(sprintf('# {\\it{}%s} cells',TYPES(k)))
        if k == 2, xlabel(xlab), else, ax(k).XTickLabel = []; end
        ax(k).PlotBoxAspectRatio = [1.3 1 1];
    end
end
