function SupPanel_A_F3(activity,...
    resp_allStim,...
    stimlabel,...
    xrange,...
    example,...
    ax)

    col = [0 0.447 0.741; 0.85 0.325 0.098];
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

    plot_vector(wt,1);
    plot_vector(mecp2,2);

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

    function plot_vector(R,line)
        % peak response of each neuron to each stimulus, smoothed across neurons
        M = smoothdata(squeeze(max(R(:,xrange,:),[],3)),1,'gaussian',5);
        axes(ax(line))
        imagesc(M)
        colormap(1-gray)
        caxis([0 0.5])
        hold on
        xline(1.5:length(xrange)-0.5,'k','LineWidth',1)
        ax(line).XAxisLocation = 'top';
        ax(line).XTick = 1:length(xrange);
        ax(line).XTickLabel = stimlabel(xrange);
        ax(line).TickLength = [0 0];
        ax(line).FontSize = 10;
        ax(line).LineWidth = 1;
        if line == 1
            ylabel({'{\bf{}wild type}','\color{black}# Neuron'},'Color',col(1,:))
        else
            ylabel({'{\bf{\it{}mecp2}-mutant}','\color{black}# Neuron'},'Color',col(2,:))
            pos = ax(line).Position;
            c = mycolorbar('\Delta{}F/F');
            c.Location = 'southoutside';
            c.Position = [pos(1) pos(2)-0.1 pos(3) 0.025];
            ax(line).Position = pos;
        end
    end
end
