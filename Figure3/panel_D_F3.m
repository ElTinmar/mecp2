function panel_D_F3(DFF_prop,Binary_prop,ax)

    col = [0 0.447 0.741; 0.85 0.325 0.098];

    [p,mx,sx,my,sy] = barplot_ranksum(ax,...
        DFF_prop{1}.allStim.prop_resp,...
        DFF_prop{2}.allStim.prop_resp,...
        '<% Resp. cells>_{fish}',...
        {'wild type','{\it{}mecp2}-mutant'},...
        col);

    disp('Proportion of responsive cells')
    disp(['WT: ' num2str(mx) ' +- ' num2str(sx)])
    disp(['MeCP2: ' num2str(my) ' +- ' num2str(sy)])
    disp(['bilateral ranksum test: ' num2str(p)])
    disp(' ')
end