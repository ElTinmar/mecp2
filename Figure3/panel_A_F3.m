function panel_A_F3(activity,...
    resp_allStim,...
    stimlabel,...
    xrange,...
    example,...
    ax)

%%% panel_A(activity,resp_allStim) 
% plot rasters of responsive cells averaged over trials, sorted 
% according to the average activity after stim onset for all stim type 
% pooled together
    
    cspan = [0 1];
    yspan = [-0.05 0.42];
    col = [0 0.447 0.741; 0.85 0.325 0.098];
    x_t = activity.wdw / activity.fps;
    ind_mecp2 = find(activity.ismecp2);
    ind_wt = find(~activity.ismecp2);

    % choose example
    ex_wt = ind_wt(example.wt);
    ex_mecp2 = ind_mecp2(example.mecp2);

    wt = get_average_resp(ex_wt);
    mecp2 = get_average_resp(ex_mecp2);
    wt_m = squeeze(nanmean(wt,1));
    wt_se = squeeze(nanstd(wt,1))./sqrt(size(wt,1));
    mecp2_m = squeeze(nanmean(mecp2,1));
    mecp2_se = squeeze(nanstd(mecp2,1))./sqrt(size(mecp2,1));

    plot_rasters(wt,1);
    plot_rasters(mecp2,2);

    function res = get_average_resp(file)
        % average over trials
        deltaFoFNormTrialAverage = squeeze(mean(activity.deltaFoFNorm{file}));
        % select responsive cells
        res = deltaFoFNormTrialAverage(:,resp_allStim{file},:);
        % reorder to have cells, stimulus, time
        res = permute(res,[2,1,3]);
        % order cells based on activity
        pooled = reshape(res,size(res,1),[]);
        % sort according to average activity
        [val,idx] = sort(mean(pooled,2),'descend');
        % remove cells that are not very active
        idx(val<prctile(val,5)) = [];
        res = res(idx,:,:);
    end

    function plot_rasters(R,line)
        for stim = 1:length(xrange)
            axes(ax(stim,line))
            imagesc(squeeze(R(:,xrange(stim),:)))
            colormap(1-gray)
            caxis(cspan)
            [~,ax(stim,line).XTick] = min(abs(x_t - [-3 0 3 6]'),[],2);
            if line == 1
                ax(stim,line).XTickLabels = [];
            else
                ax(stim,line).XTickLabels = [-3 0 3 6];
            end
            hold on
            plot([activity.onset...
                activity.onset],...
                [0 size(R,1)+1],'r:');
            if stim > 1
                ax(stim,line).YTickLabels = [];
            else
                if line==1
                    ylabel({'{\bf{}wild type}','\color{black}# Neuron'},'Color',col(1,:))
                else
                    ylabel({'{\bf{\it{}mecp2}-mutant}','\color{black}# Neuron'}, 'Color',col(2,:))
                end
            end
            if stim == length(xrange) 
                if line == 2
                    c = mycolorbar('\Delta{}F/F');
                    c.Position(1) = 0.945;
                    c.Label.Position(1) = -0.09;
                    c.Label.Position(2) = 0.5;
                    c.Label.Color = 'w';
                    c.Label.FontSize = 8;
                    xlabel('Time (s)');
                    ax(stim,line).XLabel.Position(1) = -25;
                end
            end
            ax(stim,line).LineWidth = 1;
            ax(stim,line).Units = 'centimeters'; 
            % tick length are norm relative to longest the visible axis
            norm = max(ax(stim,line).Position(3:4));
            ax(stim,line).TickLength = [0.04/norm 0.04/norm]; 
            ax(stim,line).Units = 'normalized'; 
            ax(stim,line).FontSize = 10; 
            if line==1
                title(ax(stim,line),stimlabel(xrange(stim)));
            end
        end
    end
end