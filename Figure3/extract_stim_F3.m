function res = extract_stim_F3(datafiles,time_window,fps)
%%% extract_stim_contrast(datafiles,time_window,fps)
% extract fluorescence an calcium events around each stim in a time
% window
%
% input : 
%     datafiles : data structure with the following fields
%         datafiles.celldata : contains the name of the CELLDATA files
%         datafiles.stimfile : contains the name of the STIMULATION files
%         datafiles.numfiles : number of files
%         datafiles.ismecp2  : 1 if the fish is MeCP2, 0 if WT
%     time_window : 2 elements vector (time before and after onset in sec) 
%     fps : number of frames per second for the 2-photon recordings
%
% output :
%     res.deltaFoF : cell array of DeltaF around stim for each file
%     res.deltaFoFNorm : same with average activity before onset set to 0
%     res.binary : cell array of Ca events around stim for each file
%     res.binaryNorm : same with average activity before onset set to 0
%     res.ismecp2 : 1 if the fish is MeCP2, 0 if WT
%     res.numfiles : number of files
%     res.fps : number of frames per second for the 2-photon recordings
%     res.time_window : 2 elements vector (time before and after onset in sec) 
%     res.window : time window in frames
%     res.onset : frame of stim onset during window
 
    wdw = round(time_window(1)*fps):round(time_window(2)*fps);
    onset = find(wdw==0);
    before = 1:onset-1;

    for file=1:datafiles.numfiles
        
        % load data from files
        log_progress('Calc. file %d/%d\n',file,datafiles.numfiles);
        CELL_DATA = load(datafiles.celldata{file});
        RASTER = load(datafiles.rasterdata{file});
        timings = load(datafiles.stimfile{file});
        
        deltaFoF = RASTER.deltaFoF;
        deltaFoF(logical(RASTER.movements),:) = 0;
        rasterClean = double(RASTER.raster);
        rasterClean(logical(RASTER.movements),:) = 0;
        
        n_neurons = size(deltaFoF,2);
        n_stimulus_type = 5; % 4 stim + 1 dummy stim
        n_stimulus_repetitions = 20;
        
        % use spontaneous activity before experiment 200s to 500s as dummy stim
        dummy_timing = round(linspace(200,500,n_stimulus_repetitions));

        % get stimulation onset and offset for each stimulus 
        frame_onset = zeros(n_stimulus_repetitions,n_stimulus_type);
        frame_offset = zeros(n_stimulus_repetitions,n_stimulus_type);
        for stim_type = 1:n_stimulus_type
            if stim_type == 1
                frame_onset(:,stim_type) = round(dummy_timing*fps);
                frame_offset(:,stim_type) = round((dummy_timing+1)*fps);
            else
                index = find(timings.randomOrder == (stim_type-1));
                frame_onset(:,stim_type) = ...
                    round(timings.timeStart(index)*fps);
                frame_offset(:,stim_type) = ...
                    round(timings.timeStop(index)*fps);
            end
        end
        
        % extract fluorescence and calcium events around stimulus onset
        deltaFoF_aroundStim = zeros(n_stimulus_repetitions,...
                                    n_stimulus_type,...
                                    n_neurons,...
                                    length(wdw));
                                
        CaEvents_aroundStim = zeros(n_stimulus_repetitions,...
                                    n_stimulus_type,...
                                    n_neurons,...
                                    length(wdw));
                                
        for stim_type = 1:n_stimulus_type
            for stim_rep = 1:n_stimulus_repetitions
                deltaFoF_aroundStim(stim_rep,stim_type,:,:) = ...
                    deltaFoF(frame_onset(stim_rep,stim_type)+wdw,:)';
                CaEvents_aroundStim(stim_rep,stim_type,:,:) = ...
                    rasterClean(frame_onset(stim_rep,stim_type)+wdw,:)';
            end
        end
        
        % subtract baseline before stim
        deltaFoF_aroundStimNorm = deltaFoF_aroundStim ...
            - nanmean(deltaFoF_aroundStim(:,:,:,before),4);
        CaEvents_aroundStimNorm = CaEvents_aroundStim ...
            - nanmean(CaEvents_aroundStim(:,:,:,before),4);

        res.deltaFoF{file} = deltaFoF_aroundStim;
        res.deltaFoFNorm{file} = deltaFoF_aroundStimNorm;
        res.binary{file} = CaEvents_aroundStim;
        res.binaryNorm{file} = CaEvents_aroundStimNorm;
        res.anatomy{file} = CELL_DATA.avg;
        res.cells{file} = CELL_DATA.cells;
    end
    
    res.ismecp2 = datafiles.ismecp2;
    res.numfiles = datafiles.numfiles;
    res.fps = fps;
    res.time_window = time_window;
    res.wdw = wdw;
    res.onset = onset;

     % report number of cells loaded
     ind_wt = find(~res.ismecp2);
     ind_mecp2 = find(res.ismecp2);
     cells_wt = cellfun(@(x) size(x,3), res.deltaFoF(ind_wt));
     cells_mecp2 = cellfun(@(x) size(x,3), res.deltaFoF(ind_mecp2));
     disp(['WT: ' ...
         num2str(length(cells_wt)) ...
         ' fish, ' ...
         num2str(mean(cells_wt)) ...
         ' +- ' ...
         num2str(std(cells_wt)./sqrt(length(cells_wt))) ...
         ' cells per fish'])
     disp(['MeCP2: ' ...
         num2str(length(cells_mecp2)) ...
         ' fish, ' ...
         num2str(mean(cells_mecp2)) ...
         ' +- ' ...
         num2str(std(cells_mecp2)./sqrt(length(cells_mecp2))) ...
         ' cells per fish'])
     disp(' ')
end