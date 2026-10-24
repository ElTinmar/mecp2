function [decoding,confusion,accuracy] = decoder(...
    activity,...
    time_span,...
    stim_subset,...
    nrep,...
    N_samp)
%%% 
% input
%
% output
%

    %% average activity during a time window after stim onset
    frame_span = round(activity.fps*time_span);
    after = activity.onset+1:activity.onset+frame_span;
    ind_mecp2 = find(activity.ismecp2);
    ind_wt = find(~activity.ismecp2);
    
    %% average activity after stim
    Mat = {};
    for file = 1:activity.numfiles
        % average activity over time window
        timeAvg = nanmean(activity.deltaFoF{file}(:,stim_subset,:,after),4);
        % reorganize the matrix to be n_neurons x n_stim x n_trials
        Mat{file} = permute(timeAvg,[3, 2, 1]);
    end
    
    %% decode activity
    A = {};
    C = {};
    D = {};
    for spl = 1:length(N_samp)
        log_progress('Sample size: %d\n',N_samp(spl))
        for file = 1:activity.numfiles
            for rep = 1:nrep
                % create random sample 
                ncells = size(Mat{file},1);
                samp = randsample(ncells,N_samp(spl));
                Mat_samp = Mat{file}(samp,:,:);
                [D{spl,rep,file},...
                C{spl,rep,file},...
                A{spl,rep,file}] = MLD(Mat_samp);
            end
        end
    end
    
    decoding{1} = D(:,:,ind_wt);
    decoding{2} = D(:,:,ind_mecp2);
    confusion{1} = C(:,:,ind_wt);
    confusion{2} = C(:,:,ind_mecp2);
    accuracy{1} = A(:,:,ind_wt);
    accuracy{2} = A(:,:,ind_mecp2);
end