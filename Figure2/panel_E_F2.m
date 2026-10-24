function panel_E_F2(DFF_prop,Binary_prop,xrange,ax)

    xticklab = {0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1};
    col = [0 0.447 0.741; 0.85 0.325 0.098];
    yspan = [0 2.5];
    ylab = '<Resp. duration (s)>_{fish}';
    xlab = 'Contrast';

    quantity{1} = DFF_prop{1}.perStim.fwhm;
    quantity{2} = DFF_prop{2}.perStim.fwhm;

    [p_local,p_global,anova_tbl,p_raw] = property_perStim_ranksum(...
        ax,...
        quantity,...
        yspan,...
        ylab,...
        xticklab,...
        xrange,...
        xlab,...
        col);

    disp('Response duration:')
    disp(['Local p-values, ranksum (plotted contrasts): ' num2str(p_raw(xrange))])
    disp(['    Holm-corrected: ' num2str(p_local(xrange))])
    disp('Global p-values, mixed-design ANOVA (within effects GG-corrected): ')
    disp(['    MeCP2/WT: ' num2str(p_global(1))])
    disp(['    Contrast: ' num2str(p_global(2))])
    disp(['    Interaction: ' num2str(p_global(3))])
    disp('Within-fish ANOVA table: ')
    disp(anova_tbl.Properties.Description)
    disp(anova_tbl)
end