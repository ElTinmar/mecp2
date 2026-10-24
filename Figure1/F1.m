function F1(datadir)
    rng('default')
    table_log = start_table('table_F1.txt'); %#ok<NASGU>
    
    %% load the file names ----------------------------------------------------
    datafiles = get_datafiles_F1(fullfile(datadir,'Spontaneous'));

    %% load data --------------------------------------------------------------
    fps = 1000/256;
    frames = 10 + (1:round(35*60*fps)); % remove first frames and cut at 35 min
    neuralActivity = load_data_F1(datafiles,fps,frames);

    %% compute correlations ---------------------------------------------------
    scale = 2.00429; % pixels/um
    rep_null = 500; % number of repetiotions for the null model
    corrdata = get_correlation_F1(neuralActivity,scale,rep_null);

    %% example fish -----------------------------------------------------------
    time_slice = round(15*60*fps);

    ex.wt.raster_ind = 4;
    ex.wt.corr_ind = 5;
    ex.wt.raster_frames = 2500 + (0:time_slice);
    ex.mecp2.raster_ind = 3;
    ex.mecp2.corr_ind = 2; 
    ex.mecp2.raster_frames = 200 + (0:time_slice);

    %% save results -----------------------------------------------------------
    save('fig1.mat',...
        'neuralActivity',...
        'ex',...
        'corrdata',...
        '-v7.3');

    %% plot results -----------------------------------------------------------
    f = plot_F1(...
        neuralActivity,...
        ex,...
        corrdata);

    %% save figures -----------------------------------------------------------
    save_figure(f,'fig1');

end
