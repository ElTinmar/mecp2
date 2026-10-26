function panel_A_FS1(ex,ax)
%%% one optical section: DAPI, gad1b, vglut2a and merge (vglut2a magenta, gad1b green, DAPI blue)

    if ~isfile(ex.file)
        for k = 1:4, placeholder_FS1(ax(k),{'example section not found:',ex.file}); end
        return
    end
    TOP = 99.5;  % display range: 0 to this percentile of each channel, linear
    linear = @(I) min(I / prctile(I(:),TOP),1);
    page = @(c) linear(double(imread(ex.file,3*ex.plane - 3 + c)));
    vglut = page(1);
    gad = page(2);
    dapi = page(3);

    imgs = {dapi, gad, vglut, cat(3,vglut,gad,max(vglut,dapi))};
    labels = {'DAPI','{\it{}gad1b}','{\it{}vglut2a}',''};
    label_col = {[0.3 0.5 1],[0 0.9 0],[1 0 1],'w'};
    for k = 1:4
        axes(ax(k))
        imshow(imgs{k})
        [h,w] = size(imgs{k},[1 2]);
        text(0.04*w,0.1*h,labels{k},'VerticalAlignment','baseline','Color',label_col{k},'FontSize',10)
    end
    fprintf('Panel A: %s, plane %d\n\n',ex.file,ex.plane);
end
