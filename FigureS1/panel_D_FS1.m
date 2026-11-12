function panel_D_FS1(cells,ax)
%%% dorso-ventral distribution of vglut2a (ax(1)) and gad1b (ax(2)) cells

    disp('Dorso-ventral distribution (D)')
    distribution_FS1(cells,'plane',[0 75],25,2,'Dorso-ventral axis (\mum)',ax);  % 2 um per plane
end
