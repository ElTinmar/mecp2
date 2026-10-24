function panel_B_F3(DFF_prop,Binary_prop,ax)
    
    col = [0 0.447 0.741; 0.85 0.325 0.098];

    [p,mx,sx,my,sy] = barplot_ranksum(ax,...
        nanmean(DFF_prop{1}.allStim.itg(:,2:end),2),...
        nanmean(DFF_prop{2}.allStim.itg(:,2:end),2),...
        '<{\fontsize{18}\int} \Delta{}F/F>_{fish}',...
        {'wild type','{\it{}mecp2}-mutant'},...
        col);

    % report mean, se, p-value
    disp('Integral of the response')
    disp(['WT: ' num2str(mx) ' +- ' num2str(sx)])
    disp(['MeCP2: ' num2str(my) ' +- ' num2str(sy)])
    disp(['bilateral ranksum test: ' num2str(p)])
    disp(' ')
end