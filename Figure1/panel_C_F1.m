function panel_C_F1(data,ax)

    yspan = [-0.05 0.42];
    col = [0 0.447 0.741; 0.85 0.325 0.098];

    duration{1} = get_duration(data.wt.raster, data.fps);
    duration{2} = get_duration(data.mecp2.raster, data.fps);

    [pval,m_wt,se_wt,m_mecp2,se_mecp2] = barplot_ranksum(...
        ax,...
        duration{1},...
        duration{2},...
        '<Duration (s)>_{fish}',...
        {'wild type','{\it{}mecp2}-mutant'},...
        col);

    ax.XTickLabelRotation = 0;

    % report mean, se, p-value
    disp('Duration of Ca events')
    disp(['WT: ' num2str(m_wt) ' +- ' num2str(se_wt)])
    disp(['MeCP2: ' num2str(m_mecp2) ' +- ' num2str(se_mecp2)])
    disp(['bilateral ranksum test: ' num2str(pval)])
    disp(' ')

    function duration = get_duration(rasters ,fps)
        duration = [];
        for file = 1:numel(rasters)
            X = rasters{file}';
            numcells = size(X,1);
            event_duration = cell(numcells,1);
            for c = 1:numcells
                [~,start] = find(diff([0 X(c,:) 0])==1);
                [~,stop] = find(diff([0 X(c,:) 0])==-1);
                assert(length(start)==length(stop));
                event_duration{c} = (stop-start)/fps;
            end
            duration(file)=nanmean(cellfun(@mean,event_duration));
        end
    end
end