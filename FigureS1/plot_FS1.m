function f = plot_FS1(cells,counts,ex)

    % create A4 sized figure ----------------------------------------------
    f = myfigure();

    % panel layout --------------------------------------------------------
    page_margin = 1.5;
    pos = [page_margin page_margin 21-2*page_margin 29.7-2*page_margin];
    main_panel = mypanel(f,pos);

    pA = mypanel_norm(main_panel,[0 .6 .66 .4],'a');
    pB = mypanel_norm(main_panel,[.66 .8 .34 .2],'b');
    pC = mypanel_norm(main_panel,[.66 .6 .34 .2],'c');
    pD = mypanel_norm(main_panel,[0 .2 1/3 .38],'d');
    pE = mypanel_norm(main_panel,[1/3 .2 1/3 .38],'e');
    pF = mypanel_norm(main_panel,[2/3 .2 1/3 .38],'f');

    % axes layout ---------------------------------------------------------
    side = 4.8; gap = 0.3;  % cm: square images, same horizontal and vertical spacing
    pA.Units = 'centimeters';
    origin = (pA.Position(3:4) - [0 0.6] - (2*side + gap))/2;  % centre the 2x2 grid below the 0.6 cm panel title
    pA.Units = 'normalized';
    for k = 1:4  % DAPI, gad1b (top); vglut2a, merge (bottom)
        ax_A(k) = myaxes(pA,[0 0 1 1]);
        ax_A(k).Units = 'centimeters';
        ax_A(k).InnerPosition = [origin + [(mod(k,2) == 0), (k <= 2)]*(side + gap), side, side];
    end
    ax_B = myaxes(pB,[.25 .2 .7 .7]);
    ax_C = myaxes(pC,[.25 .2 .7 .7]);
    ax_D = [myaxes(pD,[.25 .55 .7 .4]) myaxes(pD,[.25 .08 .7 .4])];
    ax_E = [myaxes(pE,[.25 .55 .7 .4]) myaxes(pE,[.25 .08 .7 .4])];
    ax_F = [myaxes(pF,[.25 .55 .7 .4]) myaxes(pF,[.25 .08 .7 .4])];

    % draw panels ---------------------------------------------------------
    panel_A_FS1(ex,ax_A);
    panel_B_FS1(counts,ax_B);
    panel_C_FS1(counts,ax_C);
    disp('Dorso-ventral distribution (D)')
    distribution_FS1(cells,'plane',[0 75],25,2,'Dorso-ventral axis (\mum)',ax_D);  % 2 um per plane
    panel_E_FS1(ax_E);
    disp('Latero-medial distribution (F), along the tectal midline')
    distribution_FS1(cells,'tectalLine',[0 1],50,1,'Latero-medial axis (norm.)',ax_F);
    col = [0 0.447 0.741; 0.85 0.325 0.098];
    text(ax_F(1),0.98,0.95,'wild type','Units','normalized','HorizontalAlignment','right',...
        'Color',col(1,:),'FontWeight','bold','FontSize',8);
    text(ax_F(1),0.98,0.84,'{\it{}mecp2}-mutant','Units','normalized','HorizontalAlignment','right',...
        'Color',col(2,:),'FontWeight','bold','FontSize',8);

    % apply axes properties ----------------------------------------------
    myaxes_properties([ax_B ax_C ax_D ax_E ax_F]);
end
