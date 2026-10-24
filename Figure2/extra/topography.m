function topography(datafiles,file,resp_perStim)
%%% topography(datafiles,resp_perStim)

    CELL_DATA = load(datafiles.celldata{file});
    cmap = jet(size(resp_perStim,2));
    try
        morpho = imadjust(CELL_DATA.avg_clean)
    catch
        morpho = imadjust(CELL_DATA.avg)
    end
    figure('units','normalized','position',[0 0 1 1]);
    for stim = 1:size(resp_perStim,2)
        subplot(4,3,stim)
        imshow(morpho);
        hold on
        for cell = resp_perStim{file,stim}
            patch(CELL_DATA.cell_per{cell}(:,1),CELL_DATA.cell_per{cell}(:,2),cmap(stim,:));
        end
        drawnow
    end
end