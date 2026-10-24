function panel_B_F1(data,ax)
    
    yspan = [-0.05 0.42];
    col = [0 0.447 0.741; 0.85 0.325 0.098];

    frequency{1} = get_frequency(data.wt.raster, data.fps);
    frequency{2} = get_frequency(data.mecp2.raster, data.fps);

    [pval,m_wt,se_wt,m_mecp2,se_mecp2] = barplot_ranksum(...
        ax,...
        frequency{1},...
        frequency{2},...
        '<Frequency (Hz)>_{fish}',...
        {'wild type','{\it{}mecp2}-mutant'},...
        col);

    ax.XTickLabelRotation = 0;

    % report mean, se, p-value
    disp('Frequency of Ca events')
    disp(['WT: ' num2str(m_wt) ' +- ' num2str(se_wt)])
    disp(['MeCP2: ' num2str(m_mecp2) ' +- ' num2str(se_mecp2)])
    disp(['bilateral ranksum test: ' num2str(pval)])
    disp(' ')

    function frequency = get_frequency(rasters ,fps)
        frequency = [];
        for file = 1:numel(rasters)
            X = rasters{file}';
            numcells = size(X,1);
            freq = cell(numcells,1);
            for c = 1:numcells
                [~,start] = find(diff([0 X(c,:) 0])==1);
                freq{c} = fps * length(start)./size(X,2);
            end
            frequency(file)=nanmean(cellfun(@mean,freq));
        end
    end 
end