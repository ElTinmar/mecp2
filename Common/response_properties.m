function [DFF, Binary] = response_properties(activity,resp_perStim,resp_allStim,time_span)
%%% [peakDeltaFoF, peakBinary] = peak_response(activity,resp_perStim,time_span)
% return peak response / integral / fwhm per stim and global
%
% input
%
% output
%
    frame_span = round(activity.fps*time_span);
    after = activity.onset:activity.onset+frame_span;
    
    ind_mecp2 = find(activity.ismecp2);
    ind_wt = find(~activity.ismecp2);

    n_stim = size(resp_perStim,2);

    DFF{1} = get_properties(ind_wt,activity.deltaFoFNorm);
    Binary{1} = get_properties(ind_wt,activity.binaryNorm);
    Binary{2} = get_properties(ind_mecp2,activity.binaryNorm);
    DFF{2} = get_properties(ind_mecp2,activity.deltaFoFNorm);

    function res = get_properties(ind,act)
        res = [];
        for i = 1:length(ind)
            % average over trials
            TrialAvg = squeeze(mean(act{ind(i)}));
            resp_cells_all = resp_allStim{ind(i)};

            for stim = 1:n_stim
                %% cells that respond to a given stim type
                resp_cells = resp_perStim{ind(i),stim};
                if length(resp_cells) == 1 % if there is only one cell
                    TrialAvgRespCellAvg = squeeze(TrialAvg(stim,resp_cells,:));
                elseif length(resp_cells) > 1
                    TrialAvgRespCellAvg = squeeze(nanmean(TrialAvg(stim,resp_cells,:)));
                else
                    TrialAvgRespCellAvg = [];
                end
                [pk,itg,fwhm] = get_prop(TrialAvgRespCellAvg,false);

                res.perStim.pk(i,stim) = pk;
                res.perStim.itg(i,stim) = itg;
                res.perStim.fwhm(i,stim) = fwhm;
                res.perStim.prop_resp(i,stim) = 100*length(resp_cells)./size(act{ind(i)},3);

                %% cells that respond for at least one stim type
                if length(resp_cells_all) == 1 % if there is only one cell
                    TrialAvgRespCellAvg = squeeze(TrialAvg(stim,resp_cells_all,:));
                elseif length(resp_cells_all) > 1
                    TrialAvgRespCellAvg = squeeze(nanmean(TrialAvg(stim,resp_cells_all,:)));
                else
                    TrialAvgRespCellAvg = [];
                end
                [pk,itg,fwhm] = get_prop(TrialAvgRespCellAvg,false);
                res.allStim.pk(i,stim) = pk;
                res.allStim.itg(i,stim) = itg;
                res.allStim.fwhm(i,stim) = fwhm;
            end
            res.allStim.prop_resp(i) = 100*length(resp_cells_all)./size(act{ind(i)},3);
        end

        function [pk,itg,fwhm] = get_prop(trace,plotflag)

            if isempty(trace)
                pk = NaN;
                itg = NaN;
                fwhm = NaN;
                return;
            end

            % peak
            if plotflag
                hf = figure;
                findpeaks(trace(after), activity.fps,...
                    'SortStr','descend',...
                    'NPeaks',1,...
                    'MinPeakDistance',5,...
                    'Annotate','extents',...
                    'MinPeakHeight',0.2,...
                    'WidthReference','halfprom');
                ylim([0 0.5])
                pause(2)
                close(hf)
            end
            
            [pk,~,fwhm] = findpeaks(trace(after), activity.fps,...
                'SortStr','descend',...
                'NPeaks',1,...
                'MinPeakDistance',5,...
                'WidthReference','halfprom');

            % integral
            itg = sum(trace(after));

            if isempty(pk)
                pk = NaN;
                fwhm = NaN;
            end
        end
    end
end