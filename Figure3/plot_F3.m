function f = plot_F3(neuralActivityAroundStim,...
    resp_allStim,...
    DFF_prop,...
    Binary_prop,...
    decoding,...
    confusion,...
    accuracy,...
    N_samp,...
    example,...
    stim_label,...
    stim_subset)

    % create A4 sized figure ----------------------------------------------
    f = myfigure();

    % panel layout --------------------------------------------------------
    page_margin = 1.5;
    pos = [page_margin page_margin 21-2*page_margin 29.7-2*page_margin];
    main_panel = mypanel(f,pos);

    panel_position = [...
        0 1/2 1 1/2; % A

        0 1/4 1/3 1/4; % B
        1/3 1/4 1/3 1/4; % C
        2/3 1/4 1/3 1/4; % D

        0 0.03 2/3 1/4;  % E
        2/3 0.03 1/3 1/4];  % F
   
    pA = mypanel_norm(main_panel,panel_position(1,:),'a');
    pB = mypanel_norm(main_panel,panel_position(2,:),'b');
    pC = mypanel_norm(main_panel,panel_position(3,:),'c');
    pD = mypanel_norm(main_panel,panel_position(4,:),'d');
    pE = mypanel_norm(main_panel,panel_position(5,:),'e');
    pF = mypanel_norm(main_panel,panel_position(6,:),'f');

    % axes layout ---------------------------------------------------------
    ind_mecp2 = find(neuralActivityAroundStim.ismecp2);
    ind_wt = find(~neuralActivityAroundStim.ismecp2);
    ncells_wt = length(resp_allStim{ind_wt(example.wt)});
    ncells_mecp2  = length(resp_allStim{ind_mecp2(example.mecp2)});
    ratio = ncells_wt/(ncells_wt+ncells_mecp2);

    w = 1/(length(stim_subset)+2);
    for i=1:length(stim_subset)
        for j=1:2
            if j==1
                ax_A(i,j) = myaxes(pA,[(w+.0425)*i-0.1 (1-ratio)+0.055  1.15*w .85*ratio]);
            elseif j==2
                ax_A(i,j) = myaxes(pA,[(w+.0425)*i-0.1 0.1 1.15*w .85*(1-ratio)]);
            end
        end
    end
    ax_B = myaxes(pB,[.25 .2 .725 .725]);
    ax_C = myaxes(pC,[.25 .2 .725 .725]);
    ax_D = myaxes(pD,[.25 .2 .725 .725]);

    ax_E(1) = myaxes(pE,[0.125 0.2 .3625 .725]);
    ax_E(2) = myaxes(pE,[.525 0.2 .3625 .725]);

    ax_F = myaxes(pF,[.25 .2 .725 .725]);

    % draw panels ---------------------------------------------------------
    panel_A_F3(neuralActivityAroundStim,...
        resp_allStim,...
        stim_label,...
        stim_subset,...
        example,...
        ax_A);
    
    panel_B_F3(DFF_prop,Binary_prop,ax_B);
    panel_C_F3(DFF_prop,Binary_prop,ax_C);
    panel_D_F3(DFF_prop,Binary_prop,ax_D);
    panel_E_F3(decoding,confusion,accuracy,stim_label,stim_subset,ax_E);
    panel_F_F3(decoding,confusion,accuracy,N_samp,ax_F);

    myaxes_properties([ax_A(:)' ax_B ax_C ax_D ax_E ax_F]);
end
