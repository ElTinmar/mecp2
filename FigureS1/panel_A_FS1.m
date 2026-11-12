function panel_A_FS1(ex,ax)
%%% one optical section: DAPI, gad1b, vglut2a and merge (vglut2a magenta, gad1b green, DAPI blue)

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
    midline(ax(1),ex,w,h);
    fprintf('Panel A: %s, plane %d\n\n',ex.file,ex.plane);
end

function midline(ax,ex,w,h)
%%% axes of panels E and F in one hemisphere.
% Red: tectal midline, from 0 (R) to 1 (C) on the axis of E, from the segmentation: each cell
% stores its closest point on the midline (VerTodo columns 11-12, as y, x).
% Yellow: towards medial (M), the end at 1 on the axis of F, drawn across the midline towards ex.corner.
    S = load(ex.segmentation,'VerTodo');
    V = S.VerTodo(S.VerTodo(:,8) == ex.seg_plane,:);
    pos = V(:,7);
    if contains(ex.segmentation,'Rigth_'), pos = 1 - pos; end  % as in load_data_FS1
    [~,k] = unique(pos);  % one point per position along the midline, in increasing order
    P = smoothdata(V(k,[12 11]),1,'gaussian',max(5,round(numel(k)/2)));

    target = ex.corner .* [w h];
    [~,i] = min(vecnorm(P - target,2,2));  % the yellow arrow crosses the midline where it is closest to target
    Y = P(i,:) + [-0.6; 0.95]*(target - P(i,:));

    hold(ax,'on')
    head = 0.025*w;
    arrow(ax,trim(P,[0.055 0.04]*w),'r',head)  % shortened, to leave room for R and C
    arrow(ax,Y,'y',head)
    medial = unit(Y(end,:) - Y(1,:));
    label(ax,P(1,:) + 0.015*w*medial,[0 0],0,'R','r')  % on the ends of the midline (0 and 1 on the axis of E), R nudged medially
    label(ax,P(end,:),[0 0],0,'C','r')
    label(ax,Y(1,:),outward(flipud(Y),head),0.015*w,'L','y')  % beyond the ends (0 and 1 on the axis of F)
    label(ax,Y(end,:),outward(Y,head),0.015*w,'M','y')
end

function Q = trim(P,g)
%%% polyline P shortened by an arc length g(1) at its start and g(2) at its end
    [s,i] = unique([0; cumsum(vecnorm(diff(P),2,2))]);
    Q = interp1(s,P(i,:),linspace(g(1),s(end) - g(2),numel(s))');
end

function u = outward(P,len)
%%% direction in which the polyline P leaves its last point, over its last len
    i = find(vecnorm(P - P(end,:),2,2) > len,1,'last');
    u = unit(P(end,:) - P(i,:));
end

function arrow(ax,P,col,len)
%%% polyline P (n x 2), ending in a filled arrowhead of length len; the line stops at its base
    tip = P(end,:);
    i = find(vecnorm(P - tip,2,2) > len,1,'last');
    d = unit(tip - P(i,:));
    n = [-d(2) d(1)];
    base = tip - len*d;
    plot(ax,[P(1:i,1); base(1)],[P(1:i,2); base(2)],'Color',col,'LineWidth',1)
    tri = tip + len*[0 0; -d + n/2; -d - n/2];
    fill(ax,tri(:,1),tri(:,2),col,'EdgeColor','none')
end

function label(ax,from,u,gap,txt,col)
%%% txt beyond the point from, in direction u, with its text box gap away from that point
    t = text(ax,from(1),from(2),txt,'Color',col,'FontSize',8,'FontWeight','bold',...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    half = t.Extent(3:4)/2;  % the box is taller than wide: its edge is farther along y than along x
    t.Position(1:2) = from + (gap + abs(u)*half')*u;
end

function u = unit(v)
    u = v/norm(v);
end
