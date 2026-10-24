function panel_F_F1(corrdata,ax)

    
    yspan = [-0.05 0.42];
    col = [0 0.447 0.741; 0.85 0.325 0.098];

    jsd_wt = get_jsd(corrdata.wt);
    jsd_mecp2 = get_jsd(corrdata.mecp2);

    [pval,m_wt,se_wt,m_mecp2,se_mecp2] = barplot_ranksum(...
        ax,...
        jsd_wt,...
        jsd_mecp2,...
        '<JSD>_{fish}',...
        {'wild type','{\it{}mecp2}-mutant'},...
        col);

    ax.XTickLabelRotation = 0;

    % report mean, se, p-value
    disp('JSD')
    disp(['WT: ' num2str(m_wt) ' +- ' num2str(se_wt)])
    disp(['MeCP2: ' num2str(m_mecp2) ' +- ' num2str(se_mecp2)])
    disp(['bilateral ranksum test: ' num2str(pval)])
    disp(' ')

    function res = get_jsd(data)
        for fish = 1:numel(data.C)
            corr = data.corr_dist(fish,:)+eps;
            corr = corr./sum(corr);

            nrep = size(data.corr_dist_null,2);
            for r = 1:nrep
                corr_null = squeeze(data.corr_dist_null(fish,r,:))+eps;
                corr_null = corr_null'./sum(corr_null);
                jsd(r,:) = sqrt(kldiv(corrdata.corr_centers,...
                    corr,...
                    corr_null,...
                    'js'));
            end

            res(fish) = mean(jsd);
        end
    end
end