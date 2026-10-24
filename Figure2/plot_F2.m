function f = plot_F2(neuralActivityAroundStim,...
    resp_allStim,...
    DFF_prop,...
    Binary_prop,...
    decoding,...
    confusion,...
    accuracy,...
    N_samp,...
    example,...
    stim_subset)

    % create A4 sized figure ----------------------------------------------
    f = myfigure();

    % panel layout --------------------------------------------------------
    page_margin = 1.5;
    pos = [page_margin page_margin 21-2*page_margin 29.7-2*page_margin];
    main_panel = mypanel(f,pos);

    panel_position = [...
        0 1/2 1 1/2; % A

        0 1/4 1/4 1/4; % B
        1/4 1/4 1/4 1/4; % C
        2/4 1/4 1/4 1/4; % D
        3/4 1/4 1/4 1/4; % E

        0 0.03 2/3 1/4;  % F
        2/3 0.03 1/3 1/4];  % G
   
    pA = mypanel_norm(main_panel,panel_position(1,:),'a');
    pB = mypanel_norm(main_panel,panel_position(2,:),'b');
    pC = mypanel_norm(main_panel,panel_position(3,:),'c');
    pD = mypanel_norm(main_panel,panel_position(4,:),'d');
    pE = mypanel_norm(main_panel,panel_position(5,:),'e');
    pF = mypanel_norm(main_panel,panel_position(6,:),'f');
    pG = mypanel_norm(main_panel,panel_position(7,:),'g');

    % axes layout ---------------------------------------------------------
    ind_mecp2 = find(neuralActivityAroundStim.ismecp2);
    ind_wt = find(~neuralActivityAroundStim.ismecp2);
    ncells_wt = length(resp_allStim{ind_wt(example.wt)});
    ncells_mecp2  = length(resp_allStim{ind_mecp2(example.mecp2)});
    ratio = ncells_wt/(ncells_wt+ncells_mecp2);

    w = 1/(length(stim_subset)+2);
    for i=1:length(stim_subset)
        for j=1:3
            if j==1
                ax_A(i,j) = myaxes(pA,[(w+0.00725)*i 0.265+(1-ratio)*0.7  0.9*w 1*ratio*0.7]);
            elseif j==2
                ax_A(i,j) = myaxes(pA,[(w+0.00725)*i 0.245 0.9*w 1*(1-ratio)*0.7]);
            elseif j==3
                ax_A(i,j) = myaxes(pA,[(w+0.00725)*i 0.06 0.9*w 1*0.2]);
            end
        end
    end

    ax_B = myaxes(pB,[.25 .25 .725 .725]);
    ax_C = myaxes(pC,[.25 .25 .725 .725]);
    ax_D = myaxes(pD,[.25 .25 .725 .725]);
    ax_E = myaxes(pE,[.25 .25 .725 .725]);

    ax_F(1) = myaxes(pF,[.25*3/8 0.25 .3625 .725]);
    ax_F(2) = myaxes(pF,[.25*3/8+0.4 0.25 .3625 .725]);
    ax_G = myaxes(pG,[.26 .25 .725 .725]);

    % draw panels ---------------------------------------------------------
    panel_A_F2(neuralActivityAroundStim,...
        resp_allStim,...
        stim_subset,...
        example,...
        ax_A);
    panel_D_F2(DFF_prop,Binary_prop,stim_subset,ax_B);
    panel_B_F2(DFF_prop,Binary_prop,stim_subset,ax_C);
    panel_C_F2(DFF_prop,Binary_prop,stim_subset,ax_D);
    panel_E_F2(DFF_prop,Binary_prop,stim_subset,ax_E);
    panel_F_F2(decoding,confusion,accuracy,stim_subset,ax_F);
    panel_G_F2(decoding,confusion,accuracy,N_samp,ax_G);

    myaxes_properties([ax_A(:)' ax_B ax_C ax_D ax_E ax_F ax_G]);
end
