function datafiles = get_datafiles_F1(directory)
%%% datafiles = get_datafiles(directory) 
% get CELLDATA files (contains the fluorescence time series, position of
% segmented neurons,...)
%
% input : string for the directory containing the data. This directory
%         should contain a MeCP2 folder and a WT folder
%
% output: data structure with the following fields
%   datafiles.celldata : contains the name of the CELLDATA files
%   datafiles.numfiles : number of files
%   datafiles.ismecp2  : 1 if the fish is MeCP2, 0 if WT

    %% get all CELLDATA files inside the directory
    celldata = dir([directory '/**/*CELL_DATA.mat']);
    datafiles.numfiles = numel(celldata);

    for i = 1:datafiles.numfiles
        datafiles.celldata{i} = fullfile(celldata(i).folder, ...
                                    celldata(i).name);
         
        % use the directory structure to find if the fish is WT or MeCP2-/-
        [~,cut] = regexp(datafiles.celldata{i},directory);
        if contains(datafiles.celldata{i}(cut+1:end),'MeCP2')
            datafiles.ismecp2(i) = 1; 
        elseif contains(datafiles.celldata{i}(cut+1:end),'WT')
            datafiles.ismecp2(i) = 0; 
        end
    end
end