function data = load_data_F1(datafiles,fps,frames)
%%%
% input
%
% output

    ind_mecp2 = find(datafiles.ismecp2);
    ind_wt = find(~datafiles.ismecp2);

    log_progress('Loading wt\n')
    data.wt = get_data(ind_wt);
    log_progress('Loading mecp2\n')
    data.mecp2 = get_data(ind_mecp2);
    data.fps = fps;

    % report number of cells loaded
    cells_wt = cellfun(@(x) size(x,2), data.wt.deltaFoF);
    cells_mecp2 = cellfun(@(x) size(x,2), data.mecp2.deltaFoF);
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

    function res = get_data(ind)
        for num = 1:length(ind)
            log_progress('%d/%d\n',num,length(ind))
            CELL_DATA = load(datafiles.celldata{ind(num)});
            deltaFoF = CELL_DATA.deltaFoF;
            raster = single(CELL_DATA.raster);
            raster(logical(CELL_DATA.movements),:) = 0;
            deltaFoF(logical(CELL_DATA.movements),:) = 0;
            orig_position = CELL_DATA.orig_position;
            deltaFoF_ori = deltaFoF;
            deltaFoF(~logical(raster)) = 0;
    
            % remove silent cells
            remove_cells = unique(find(sum(raster,1)==0));
            raster(:,remove_cells) = [];
            deltaFoF(:,remove_cells) = [];
            deltaFoF_ori(:,remove_cells) = [];
            orig_position(remove_cells,:) = [];
    
            res.deltaFoF_ori{num} = deltaFoF_ori(frames,:);
            res.deltaFoF{num} = deltaFoF(frames,:);
            res.raster{num} = raster(frames,:);
            res.orig_position{num} = orig_position;
            res.movements{num} = CELL_DATA.movements(frames);
        end
    end
end