function panel_G_F2(decoding,confusion,accuracy,N_samp,ax)

    n_stim = size(decoding{1}{1},1);
    n_trials = size(decoding{1}{1},2);
    totstim = n_stim * n_trials;
    col = [0 0.447 0.741; 0.85 0.325 0.098];
    ylo = 100/n_stim-2;
    yhi = 35;
    
	% chance level (binomial distribution)
    upper_bound = 100*find(binocdf(1:totstim,totstim,1/n_stim,'upper')<0.05,1,'first')/totstim; 
    disp('Decoder:')
    disp(['Chance level (95th percentile of the binomial distribution): ', num2str(upper_bound), '%']);
    
    % stats
    % average repetitions within each fish: n_fish x n_samp, fish are the unit
    accuracy_wt = squeeze(mean(cell2mat(accuracy{1}),2))';
    accuracy_mecp2 = squeeze(mean(cell2mat(accuracy{2}),2))';
    pVal = arrayfun(@(s) ranksum(accuracy_wt(:,s),accuracy_mecp2(:,s)), 1:length(N_samp));
    pHolm = holm(pVal);
    [p_anova,anova_tbl] = mixed_anova(accuracy_wt,accuracy_mecp2);

    disp('Mixed-design ANOVA, genotype x ensemble size (within effects GG-corrected):')
    disp(['    MeCP2/WT: ' num2str(p_anova(1))])
    disp(['    Ensemble size: ' num2str(p_anova(2))])
    disp(['    Interaction: ' num2str(p_anova(3))])
    disp(anova_tbl.Properties.Description)
    disp(anova_tbl)
    sem = @(a) std(a,[],1,'omitnan')./sqrt(size(a,1));
    disp(['Accuracy (%), mean +- SEM across fish, ensemble sizes ' mat2str(N_samp) ':'])
    fprintf('    WT:    %s\n', sprintf('%.2f +- %.2f   ', 100*[mean(accuracy_wt,1,'omitnan'); sem(accuracy_wt)]))
    fprintf('    MeCP2: %s\n', sprintf('%.2f +- %.2f   ', 100*[mean(accuracy_mecp2,1,'omitnan'); sem(accuracy_mecp2)]))
    disp(' ')
    disp('Post-hoc: bilateral ranksum test per ensemble size, across fish (repetitions averaged per fish):')
    disp('p-values: ')
    disp(pVal)
    disp('Holm-corrected: ')
    disp(pHolm)

    % shade the range of sample sizes from which all larger ones are significant
    start = find(pHolm>=0.05,1,'last') + 1;
    if isempty(start), start = 1; end

    axes(ax)
    if start <= length(N_samp)
        patch([N_samp(start) N_samp(end) N_samp(end) N_samp(start) N_samp(start)],...
            [ylo ylo yhi yhi ylo],...
            [0.3 0.3 0.3],...
            'FaceAlpha',0.3,...
            'EdgeColor','none')
    end
    hold on
    errorbar(N_samp,...
        100*nanmean(accuracy_wt,1),...
        100*nanstd(accuracy_wt,[],1)./sqrt(size(accuracy_wt,1)),...
        'Color',col(1,:),...
        'LineWidth',2)
    errorbar(N_samp,...
        100*nanmean(accuracy_mecp2,1),...
        100*nanstd(accuracy_mecp2,[],1)./sqrt(size(accuracy_mecp2,1)),...
        'Color',col(2,:),...
        'LineWidth',2)

    % legend
    text(15, 31, 'wild type', 'FontSize', 8, 'FontWeight', 'bold', 'Color', col(1,:));
    text(15, 33, '{\it{}mecp2}-mutant', 'FontSize', 8, 'FontWeight', 'bold', 'Color', col(2,:));
    
    plot([0 N_samp(end)],[100/n_stim 100/n_stim],'k','LineWidth',2);
    plot([0 N_samp(end)],[upper_bound upper_bound],':k','LineWidth',1)
    ylabel('<Decoding performance (%)>_{fish}')
    xlabel('Number of neurons')
    ylim([ylo yhi])
    xlim([-10 515])
    ax.XTick = N_samp(gradient(N_samp)>N_samp(end)/5);
end