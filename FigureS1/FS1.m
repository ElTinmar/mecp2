function FS1(datadir)
    rng('default')
    table_log = start_table('table_S1.txt'); %#ok<NASGU>

    %% load data --------------------------------------------------------------
    cells = load_data_FS1(fullfile(datadir,'Immunostaining','ALL.mat'));
    counts = count_cells_FS1(cells);

    %% example optical section (panel A) ------------------------------------
    ex.file = fullfile(datadir,'Immunostaining','WT','Fish 7-wt-6dpf_060318.tif');
    ex.plane = 33;  % z in the stack, 1-based (pages 3*plane-2:3*plane = vglut2a, gad1b, DAPI)

    %% plot results and statistics -------------------------------------------
    f = plot_FS1(cells,counts,ex);

    %% save results -----------------------------------------------------------
    save('figS1.mat','counts','ex');
    save_figure(f,'figS1');
end
