function datafiles = get_datafiles_F3(directory)
%%% datafiles = get_datafiles(directory) 
% get CELLDATA files (contains the fluorescence time series, position of
% segmented neurons,...) and STIMULATION file (contains the timings and
% nature of stimuli)
%
% input : string for the directory containing the data. This directory
%         should contain a MeCP2 folder and a WT folder
%
% output: data structure with the following fields
%   datafiles.celldata : contains the name of the CELLDATA files
% 	datafiles.stimfile : contains the name of the STIMULATION files
%   datafiles.numfiles : number of files
%   datafiles.ismecp2  : 1 if the fish is MeCP2, 0 if WT

    %% get all MAT files inside the directory
    allmat = dir([directory '/**/*.mat']);
    
    %% use regular expression to find STIMULATION file and CELLDATA files
    re_stim = '[0-9]{2}_1-95Hz.mat'; 
    position_match = regexp({allmat.name},re_stim);
    index_match = ~cellfun('isempty', position_match);
    stimfile = allmat(index_match);

    re_celldata = 'ALL_CELLS.mat'; 
    position_match = regexp({allmat.name},re_celldata);
    index_match = ~cellfun('isempty', position_match);
    celldata = allmat(index_match);

    re_rasterdata = 'RASTER.mat'; 
    position_match = regexp({allmat.name},re_rasterdata);
    index_match = ~cellfun('isempty', position_match);
    rasterdata = allmat(index_match);

    %% verify that we have the same number of CELLDATA and STIMULATION file
    assert(numel(stimfile) == numel(celldata))
    datafiles.numfiles = numel(stimfile);

    for i = 1:datafiles.numfiles
        datafiles.celldata{i} = fullfile(celldata(i).folder, ...
                                    celldata(i).name);
        datafiles.rasterdata{i} = fullfile(rasterdata(i).folder, ...
                                    rasterdata(i).name);
        datafiles.stimfile{i} = fullfile(stimfile(i).folder, ...
                                    stimfile(i).name);
         
        % use the directory structure to find if the fish is WT or MeCP2-/-
        [~,cut] = regexp(datafiles.celldata{i},directory);
        if contains(datafiles.celldata{i}(cut+1:end),'MeCP2')
            datafiles.ismecp2(i) = 1; 
        elseif contains(datafiles.celldata{i}(cut+1:end),'WT')
            datafiles.ismecp2(i) = 0; 
        end
    end
end