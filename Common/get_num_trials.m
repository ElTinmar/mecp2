function num_trial = get_num_trials(boot,prop_resp,n_trials)

    % use dummy stim 1 to assess the fraction of false positives
    dummy = 1;
    for file = 1:boot.numfiles
        n_cells = size(boot.respTrials{file},2);
        for t = 1:n_trials
            resp{file,t} = find(boot.respTrials{file}(dummy,:) >= t);
            prop(file,t) = length(resp{file,t})/n_cells;
        end 
    end

    ind_wt = find(~boot.ismecp2);
    ind_mecp2 = find(boot.ismecp2);
    prop_avg = max(nanmean(prop(ind_mecp2,:)), nanmean(prop(ind_wt,:)));
    num_trial = find(prop_avg <= prop_resp,1,'first');
end