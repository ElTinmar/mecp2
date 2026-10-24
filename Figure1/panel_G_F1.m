function panel_G_F1(corrdata,ax)
    
    yspan = [-0.005 0.06];
    xspan = [0 240];
    col = [0 0.447 0.741; 0.85 0.325 0.098];

    wt_m = nanmean(corrdata.wt.corr_distance);
    wt_se = nanstd(corrdata.wt.corr_distance)./sqrt(size(corrdata.wt.corr_distance,1));
    mecp2_m = nanmean(corrdata.mecp2.corr_distance);
    mecp2_se = nanstd(corrdata.mecp2.corr_distance)./sqrt(size(corrdata.mecp2.corr_distance,1));
    
    axes(ax)
    hold on

    pVal = [];
    for d = 1:length(corrdata.dist_centers)
        pVal(d) = ranksum(corrdata.wt.corr_distance(:,d),...
                    corrdata.mecp2.corr_distance(:,d));
    end
    pHolm = holm(pVal);
    disp('Correlation vs distance, bilateral ranksum test per distance bin:')
    disp(pVal)
    disp('Holm-corrected: ')
    disp(pHolm)

    p = plot(corrdata.dist_centers,wt_m);
    p.LineWidth = 1;
    p.Color = col(1,:);
    xx = [corrdata.dist_centers fliplr(corrdata.dist_centers)];
    yy = [wt_m+wt_se fliplr(wt_m-wt_se)]; 
    h = patch(xx,yy,col(1,:));
    h.EdgeColor = 'none';
    h.FaceAlpha = 0.3;

    p1 = plot(corrdata.dist_centers,mecp2_m);
    p1.LineWidth = 1;
    p1.Color = col(2,:);
    xx = [corrdata.dist_centers fliplr(corrdata.dist_centers)];
    yy = [mecp2_m+mecp2_se fliplr(mecp2_m-mecp2_se)]; 
    h = patch(xx,yy,col(2,:));
    h.EdgeColor = 'none';
    h.FaceAlpha = 0.3;

    % add p-values to the plot
    signif = find(pHolm<=0.05);
    for  i = signif
        X = corrdata.dist_centers(i) + [-15 15 15 -15];
        Y = [yspan(1) yspan(1) yspan(2) yspan(2)];
        h = fill(X,Y,[0.3 0.3 0.3],'EdgeColor','none');
        h.FaceAlpha = 0.3;
    end

    %lgd = legend([p p1],{'wild type','{\it{}mecp2}-mutant'},'FontSize',8,'Box','off');
    %lgd.ItemTokenSize = [5 5];
    %lgd.Position(1) = 0.4;
    text(230, 0.0455, 'wild type', 'FontSize', 8, 'FontWeight', 'bold', 'Color', col(1,:), 'HorizontalAlignment', 'right');
    text(230, 0.0525, '{\it{}mecp2}-mutant', 'FontSize', 8, 'FontWeight', 'bold', 'Color', col(2,:), 'HorizontalAlignment', 'right');

    xlim(xspan);
    ylim(yspan);
    ax.XTick = [30 90 150 210];
    ylabel('<Correlation>_{fish}')
    mu = char(hex2dec('0x03BC'));
    xlabel("Distance ("+mu+"m)")
end