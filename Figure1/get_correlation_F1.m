function corrdata = get_correlation_F1(activity,scale,rep_null)
    
    % create bin edges and centers to estimate the distribution
    corr_edges = linspace(-1,1,1001);
    d = diff(corr_edges)/2;
    corr_centers = corr_edges(1:end-1)+d;
    corrdata.corr_edges = corr_edges;
    corrdata.corr_centers = corr_centers;

    dist_edges = 0:30:240;
    d = diff(dist_edges)/2;
    dist_centers = dist_edges(1:end-1)+d;
    corrdata.dist_edges = dist_edges;
    corrdata.dist_centers = dist_centers;

    log_progress('wild type\n')
    corrdata.wt = compute_corr(activity.wt);
    log_progress('MeCP2\n')
    corrdata.mecp2 = compute_corr(activity.mecp2);

    function res = compute_corr(data)
        for file = 1:numel(data.deltaFoF)
            log_progress('%d/%d\n',file,numel(data.deltaFoF))

            % remove movements artifacts to compute the correlations
            cutFrames = find(data.movements{file});
            deltaFoFcut = data.deltaFoF{file};
            deltaFoFcut(cutFrames,:) = [];

            % compute correlations
            C = corr(deltaFoFcut,'type','Pearson');
            upper_tri = triu(ones(size(C),'logical'),1);
            ind = upper_tri(:);
            res.C{file} = C;
            res.corr_dist(file,:) = histcounts(...
                C(ind),corr_edges,...
                'Normalization','probability');

            % compute distance between cells
            res.D{file} = squareform(pdist(data.orig_position{file}));
            res.D{file} = res.D{file}./scale; % convert to microns

            % put NaNs on the diagonal
            ind_diag = find(eye(size(res.C{file}),'logical'));
            res.C{file}(ind_diag) = NaN;
            res.D{file}(ind_diag) = NaN;

            % correlation vs distance 
            dbin = discretize(res.D{file}(ind),dist_edges);
            res.corr_distance(file,:) = splitapply(@(x) nanmean(x),C(ind),dbin);

            % compute null-model correlations
            log_progress('Computing null-model\n')
            n_samples = size(deltaFoFcut,1);
            n_neurons = size(deltaFoFcut,2);
            for r = 1:rep_null
                % print progress 
                if mod(r,rep_null/10)==0
                    log_progress('%g%% ',100*r/rep_null)
                end

                % compute correlations null model
                if n_samples>=n_neurons
                    % do not allow replacement
                    shifts = randsample(n_samples,n_neurons,false); 
                else
                    warning('n_samples < n_neurons, sampling with replacement')
                    shifts = randsample(n_samples,n_neurons,true);
                end
                deltaFoF_shifted = zeros(size(deltaFoFcut));
                for j=1:n_neurons
                    deltaFoF_shifted(:,j) = circshift(deltaFoFcut(:,j),shifts(j));
                end
                C_null = corr(deltaFoF_shifted,'type','Pearson');
                res.corr_dist_null(file,r,:) = histcounts(...
                C_null(ind),corr_edges,...
                'Normalization','probability');
            end 
        end
    end
end