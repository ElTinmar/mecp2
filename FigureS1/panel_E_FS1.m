function panel_E_FS1(cells,ax)
%%% rostro-caudal distribution: position along the tectal midline

    disp('Rostro-caudal distribution (E), along the tectal midline')
    distribution_FS1(cells,'rostroCaudal',[0 1],50,1,'Rostro-caudal axis (norm.)',ax);
    for a = ax, ylabel(a,''), end  % y-labels shown in D only
end
