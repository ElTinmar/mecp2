function [resp_allStim, resp_perStim] = get_resp_cells(boot,thresh_num_trial)
%%% [resp_allStim, resp_perStim] = get_resp_cells(boot,thresh_num_trial)
% Find cells responsive to the stimulus for each stim type (resp_perStim)
% or for at least one stim type (resp_allStim)
%
% input
%     boot : data structure with fields
%         boot.stat : bootstraped difference of means
%         boot.prc : 5th percentile of the difference of means
%         boot.resp : for each cell/stim/type/trial 1 if there is a response 
%         boot.respTrials : number of responsive trials for each cell/stim type
%         boot.numfiles : number of files
%         boot.ismecp2 : 1 if the fish is MeCP2, 0 if WT
%     thresh_num_trial : number of responsive trials for a cell to be deemed
%       responsive for a given stim type
% output
%     resp_allStim : cells responsive to at least <thresh_num_trial> trials of one stim
%       type
%     resp_perStim : cells responsive to at least <thresh_num_trial> trials for each stim
%       type

    for file = 1:boot.numfiles
        n_stimulus_type = size(boot.respTrials{file},1);
        for stim_type = 1:n_stimulus_type
            resp_perStim{file,stim_type} = ...
                find(boot.respTrials{file}(stim_type,:) >= thresh_num_trial);
        end
        % stim 1 is a dummy stimulation don't use it
        resp_allStim{file} = unique(horzcat(resp_perStim{file,2:end}));
    end

    % report number of cells 
    ind_wt = find(~boot.ismecp2);
    ind_mecp2 = find(boot.ismecp2);
    cells_wt = cellfun(@(x) size(x,2), resp_allStim(ind_wt));
    cells_mecp2 = cellfun(@(x) size(x,2), resp_allStim(ind_mecp2));
    disp('# of responsive cells at least one stim type: ')
    fprintf(['WT: ------------------------------------\n' ...
        num2str(length(cells_wt)) ...
        ' fish, ' ...
        num2str(mean(cells_wt)) ...
        ' +- ' ...
        num2str(std(cells_wt)./sqrt(length(cells_wt))) ...
        ' cells per fish \nTotal number of cells: ' ...
        num2str(sum(cells_wt))])
    fprintf(['\nMeCP2: ---------------------------------\n' ...
        num2str(length(cells_mecp2)) ...
        ' fish, ' ...
        num2str(mean(cells_mecp2)) ...
        ' +- ' ...
        num2str(std(cells_mecp2)./sqrt(length(cells_mecp2))) ...
        ' cells per fish \nTotal number of cells: ' ...
        num2str(sum(cells_mecp2))])
    fprintf('\n\n')
end