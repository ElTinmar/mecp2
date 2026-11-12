function panel_F_FS1(cells,ax)
%%% latero-medial distribution: distance to the tectal midline, with the genotype legend

    disp('Latero-medial distribution (F), distance to the tectal midline')
    distribution_FS1(cells,'lateroMedial',[0 1],36,1,'Latero-medial axis (norm.)',ax);
    for a = ax, ylabel(a,''), end  % y-labels shown in D only
    col = [0 0.447 0.741; 0.85 0.325 0.098];
    text(ax(1),0.97,0.93,'wild type','Units','normalized','HorizontalAlignment','right',...
        'Color',col(1,:),'FontWeight','bold','FontSize',8);
    text(ax(1),0.97,0.82,'{\it{}mecp2}-mutant','Units','normalized','HorizontalAlignment','right',...
        'Color',col(2,:),'FontWeight','bold','FontSize',8);
end
