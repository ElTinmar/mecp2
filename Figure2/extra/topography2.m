function topography2(datafiles,file,resp_allStim,boot)
%%% topography(datafiles,resp_perStim)

    CELL_DATA = load(datafiles.celldata{file});
    cmap = jet(size(boot.respTrials{file},1));
    figure('units','normalized','position',[0 0 1 1]);
    
    imshow(imadjust(CELL_DATA.avg_clean));
    hold on
    for cell = resp_allStim{file}
        [~,az] = max(boot.respTrials{file}(:,cell));
        patch(CELL_DATA.cell_per{cell}(:,1),CELL_DATA.cell_per{cell}(:,2),cmap(az,:));
    end
    drawnow
end