function f = plot_SupF3(...
    activity,...
    resp_allStim,...
    stimlabel,...
    stim_subset,...
    stim_azimuth,...
    example,...
    confusion,...
    decoding,...
    N_samp)

    % create A4 sized figure ----------------------------------------------
    f = myfigure();

    % panel layout --------------------------------------------------------
    page_margin = 1.5;
    pos = [page_margin page_margin 21-2*page_margin 29.7-2*page_margin];
    main_panel = mypanel(f,pos);

    panel_position = [...
        0 2/3 1/2 1/3;
        1/2+0.05 3/4 1/2 1/4
    ]; 
   
    pA = mypanel_norm(main_panel,panel_position(1,:),'A');
    pB = mypanel_norm(main_panel,panel_position(2,:),'B');

    % axes layout ---------------------------------------------------------
    ind_mecp2 = find(activity.ismecp2);
    ind_wt = find(~activity.ismecp2);
    ncells_wt = length(resp_allStim{ind_wt(example.wt)});
    ncells_mecp2  = length(resp_allStim{ind_mecp2(example.mecp2)});
    h_wt = ncells_wt/max(ncells_wt,ncells_mecp2);
    h_mecp2 = ncells_mecp2/max(ncells_wt,ncells_mecp2);

    h = 0.7*[h_wt h_mecp2];
    for j = 1:2
        ax_A(j) = myaxes(pA,[0.22+0.46*(j-1) 0.88-h(j) 0.28 h(j)]);
        ax_A(j).PlotBoxAspectRatioMode = 'auto';
    end
    ax_B = myaxes(pB,[.1 .1 .85 .85]);

    % draw panels ---------------------------------------------------------
    SupPanel_A_F3(activity,...
        resp_allStim,...
        stimlabel,...
        stim_subset,...
        example,...
        ax_A);
    SupPanel_B_F3(N_samp,confusion,decoding,stim_azimuth,ax_B);

    myaxes_properties([ax_A(:)' ax_B]);
end
