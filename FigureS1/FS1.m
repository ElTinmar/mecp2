function FS1(datadir)
    rng('default')
    table_log = start_table('table_S1.txt'); %#ok<NASGU>

    %% load data --------------------------------------------------------------
    cells = load_data_FS1(fullfile(datadir,'Immunostaining','Normalized_vGlu1.1-Gad0.7-DAPI1.1'));
    counts = count_cells_FS1(cells);

    %% example optical section (panel A) ------------------------------------
    ex.file = fullfile(datadir,'Immunostaining','WT','Fish 7-wt-6dpf_060318.tif');
    ex.plane = 33;  % z in the stack, 1-based (pages 3*plane-2:3*plane = vglut2a, gad1b, DAPI)
    ex.segmentation = fullfile(datadir,'Immunostaining','Normalized_vGlu1.1-Gad0.7-DAPI1.1',...
        'Rigth_VerTodo_Fish 7-wt-6dpf_060318.mat');  % hemisphere whose midline is drawn
    ex.seg_plane = ex.plane - 1;  % plane in the segmentation: z - Tectum_Begin (1 for this larva)
    ex.corner = [0.482 0.654];  % intersection of the commissures, where both tecta meet the cerebellum, fraction of the image width and height

    %% plot results and statistics -------------------------------------------
    f = plot_FS1(cells,counts,ex);

    %% save results -----------------------------------------------------------
    save('figS1.mat','counts','ex');
    save_figure(f,'figS1');
end
