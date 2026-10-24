function f = plot_F1(neuralActivity,ex,corrdata,rep_null)
    %TODO tick length 

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
        2/4 1/4 2/4 1/4; % D
        0 0.05 2/4 1/4;  % E
        2/4 .05 1/4 1/4; % F
        3/4 .05 1/4 1/4];  % G
   
    pA = mypanel_norm(main_panel,panel_position(1,:),'a');
    pB = mypanel_norm(main_panel,panel_position(2,:),'b');
    pC = mypanel_norm(main_panel,panel_position(3,:),'c');
    pD = mypanel_norm(main_panel,panel_position(4,:),'d');
    pE = mypanel_norm(main_panel,panel_position(5,:),'e');
    pF = mypanel_norm(main_panel,panel_position(6,:),'f');
    pG = mypanel_norm(main_panel,panel_position(7,:),'g');

    % axes layout ---------------------------------------------------------
    ax_A(1) = myaxes(pA,[0.1 7/10 0.85 3/10]);
    ax_A(2) = myaxes(pA,[0.1 5.75/10 0.85 1/10]);
    ax_A(3) = myaxes(pA,[0.1 2.25/10 0.85 3/10]);
    ax_A(4) = myaxes(pA,[0.1 1/10 0.85 1/10]);
    
    ax_B = myaxes(pB,[.25 .25 .725 .725]);
    ax_C = myaxes(pC,[.25 .25 .725 .725]);
    
    ax_D(1) = myaxes(pD,[.125 .25 .3625 .725]);
    ax_D(2) = myaxes(pD,[.525 .25 .3625 .725]);
    ax_E(1) = myaxes(pE,[.125 .25 .3625 .725]);
    ax_E(2) = myaxes(pE,[.525 .25 .3625 .725]);
    
    ax_F = myaxes(pF,[.25 .25 .725 .725]);
    ax_G = myaxes(pG,[.25 .25 .725 .725]);

    % draw panels ---------------------------------------------------------
    panel_A_F1(neuralActivity,ex,ax_A);
    panel_B_F1(neuralActivity,ax_B);
    panel_C_F1(neuralActivity,ax_C);
    panel_D_F1(neuralActivity,ex,ax_D);
    panel_E_F1(neuralActivity,corrdata,ax_E);
    panel_F_F1(corrdata,ax_F);
    panel_G_F1(corrdata,ax_G);

    % apply axes properties ----------------------------------------------
    myaxes_properties([ax_A ax_B ax_C ax_D ax_E ax_F ax_G]);
end