function f = plot_FS1(cells,counts,ex)

    % create A4 sized figure ----------------------------------------------
    f = myfigure();

    % panel layout --------------------------------------------------------
    page_margin = 1.5;
    pos = [page_margin page_margin 21-2*page_margin 29.7-2*page_margin];
    main_panel = mypanel(f,pos);

    pA = mypanel_norm(main_panel,[0 .55 2/3 .45],'a');
    pB = mypanel_norm(main_panel,[2/3 .775 1/3 .225],'b');
    pC = mypanel_norm(main_panel,[2/3 .55 1/3 .225],'c');
    pD = mypanel_norm(main_panel,[0 .16 1/3 .38],'d');
    pE = mypanel_norm(main_panel,[1/3 .16 1/3 .38],'e');
    pF = mypanel_norm(main_panel,[2/3 .16 1/3 .38],'f');

    % axes layout ---------------------------------------------------------
    x = [.17 .78];  % B-F: left margin and width in their panel (the plot box fills the width)
    gap = 0.3;  % cm, between images and before panel B
    set([pA pB],'Units','centimeters');
    left = x(1)*pos(3)/3;  % left-aligned with the plots of D
    side = (pA.Position(3) - left - 2*gap)/2;  % square images filling the width of A
    top = (1 - .10 - .80)*pB.Position(4);  % top-aligned with the plot of B (axes of B below)
    origin = [left, pA.Position(4) - top - (2*side + gap)];
    set([pA pB],'Units','normalized');
    for k = 1:4  % DAPI, gad1b (top); vglut2a, merge (bottom)
        ax_A(k) = myaxes(pA,[0 0 1 1]);
        ax_A(k).Units = 'centimeters';
        ax_A(k).InnerPosition = [origin + [(mod(k,2) == 0), (k <= 2)]*(side + gap), side, side];
    end
    ax_B = myaxes(pB,[x(1) .10 x(2) .80]);
    ax_C = myaxes(pC,[x(1) .10 x(2) .80]);
    ax_D = [myaxes(pD,[x(1) .53 x(2) .44]) myaxes(pD,[x(1) .06 x(2) .44])];
    ax_E = [myaxes(pE,[x(1) .53 x(2) .44]) myaxes(pE,[x(1) .06 x(2) .44])];
    ax_F = [myaxes(pF,[x(1) .53 x(2) .44]) myaxes(pF,[x(1) .06 x(2) .44])];

    % draw panels ---------------------------------------------------------
    panel_A_FS1(ex,ax_A);
    panel_B_FS1(counts,ax_B);
    panel_C_FS1(counts,ax_C);
    panel_D_FS1(cells,ax_D);
    panel_E_FS1(cells,ax_E);
    panel_F_FS1(cells,ax_F);

    % B, C: square plot boxes, right edge aligned with the plot box of F ---
    set([ax_B ax_C ax_F(1)],'Units','centimeters');  % B, C and F panels share the same x origin
    pf = ax_F(1).InnerPosition;
    ar = ax_F(1).PlotBoxAspectRatio;
    right = pf(1) + (pf(3) + min(pf(3),pf(4)*ar(1)/ar(2)))/2;
    for a = [ax_B ax_C]
        p = a.InnerPosition;
        s = min(p(3:4));
        a.InnerPosition = [right - s, p(2) + (p(4) - s)/2, s, s];
    end
    set([ax_B ax_C ax_F(1)],'Units','normalized');

    % apply axes properties ----------------------------------------------
    myaxes_properties([ax_B ax_C ax_D ax_E ax_F]);
end
