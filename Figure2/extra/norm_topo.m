function norm_topo(activity,datafiles,resp_perStim)
%%% topography(datafiles,resp_perStim)

    ind_mecp2 = find(activity.ismecp2);
    ind_wt = find(~activity.ismecp2);
    cmap = jet(size(resp_perStim,2));
    %gather_results(ind_mecp2)
    gather_results(ind_wt)
    
    function gather_results(ind)
        f = figure;
        for file = ind
            CELL_DATA = load(datafiles.celldata{file});
            for stim = 1:size(resp_perStim,2)
                for cell = resp_perStim{file,stim}
                    scatter(CELL_DATA.orig_position(cell,1)+10*rand(),CELL_DATA.orig_position(cell,2)+10*rand(),5,cmap(stim,:),'filled');
                    hold on
                end
            end
        end
        ax = gca();
        ax.Color = 'k';
    end
end