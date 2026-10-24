function res = bootstrap_diffBeforeAfterStim(activity,time_span,method,thresh)
%%% res = bootstrap_diffBeforeAfterStim(activity,time_span)
% bootstrap the difference of activity between time_span seconds before the 
% onset and time_span seconds after the onset of stimulation (sample
% time points with replacement). A cell is deemed responsive to a trial 
% if the 5th prctile of the bootstraped distribution for that trial 
% is positive.
%
% input
%     activity : data structure with the following fields
%         activity.deltaFoF : cell array of DeltaF around stim for each file
%         activity.deltaFoFNorm : same with average activity before onset set to 0
%         activity.binary : cell array of Ca events around stim for each file
%         activity.binaryNorm : same with average activity before onset set to 0
%         activity.ismecp2 : 1 if the fish is MeCP2, 0 if WT
%         activity.numfiles : number of files
%         activity.fps : number of frames per second for the 2-photon recordings
%         activity.time_window : 2 elements vector (time before and after onset in sec) 
%         activity.window : time window in frames
%         activity.onset : frame of stim onset during window
%     time_span : duration of the time window before and after stim for the
%       bootstrap
%
% output
%     res.stat : bootstraped difference of means
%     res.prc : 5th percentile of the difference of means
%     res.resp : for each cell/stim/type/trial 1 if there is a response 
%     res.respTrials : number of responsive trials for each cell/stim type
%     res.numfiles : number of files
%     res.ismecp2 : 1 if the fish is MeCP2, 0 if WT

    
    frame_span = round(activity.fps*time_span);
    before = activity.onset-frame_span:activity.onset-1;
    after = activity.onset+1:activity.onset+frame_span;

    for file=1:activity.numfiles
        log_progress('Calc. file %d/%d\n',file,activity.numfiles);
        if strcmp(method,'binary')
            X = activity.binary{file}(:,:,:,before); 
            Y = activity.binary{file}(:,:,:,after); 
        else
            X = activity.deltaFoF{file}(:,:,:,before); 
            Y = activity.deltaFoF{file}(:,:,:,after); 
        end
        
        % bootstrp samples from the rows with replacement -> permute 
        % so that rows represent time, and average over time for bootfun
        X_t = permute(X,[4 1 2 3]);
        Y_t = permute(Y,[4 1 2 3]);
        B = bootstrp(1000,@(x,y) squeeze(mean(y)-mean(x)),X_t,Y_t);
        
        % reshape into a matrix
        res.stat{file} = reshape(B,1000,size(X,1),size(X,2),size(X,3));
        res.prc{file} = squeeze(prctile(res.stat{file},thresh));
        res.resp{file} = res.prc{file} > 0;
        res.respTrials{file} = squeeze(sum(res.resp{file}));
    end
    res.numfiles = activity.numfiles;
    res.ismecp2 = activity.ismecp2;
end