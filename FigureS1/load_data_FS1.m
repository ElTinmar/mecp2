function cells = load_data_FS1(datadir)
%%% cells = load_data_FS1(datadir)
% one row per segmented cell, from the <side>_VerTodo_Fish <n>-<genotype>-*.mat files of Segment_Immuno.m.
% VerTodo columns used: 7 position along the tectal midline (0-1, rostro-caudal), 8 plane,
% 13 signed distance to the midline (px, latero-medial), 14 cell type (1 vglut2a, 2 gad1b, 3 DAPI only).

    PLANE_MAX = 86;                 % last imaging plane
    DIST_MAX = 110;                 % px, |distance to the midline| mapped onto the latero-medial axis [0 1]
    FLIPPED = ["wt1","wt2","wt3"];  % larvae imaged bottom-up
    TYPES = ["vglut2a","gad1b","dapi"];

    files = dir(fullfile(datadir,'*_VerTodo_*.mat'));
    assert(~isempty(files),'no VerTodo files in %s',datadir);
    parts = cell(numel(files),1);
    for k = 1:numel(files)
        tok = regexp(files(k).name,'^(Left|Rigth)_VerTodo_Fish (\d+)-(\w+)-','tokens','once');
        assert(numel(tok) == 3,'unexpected file name: %s',files(k).name);
        S = load(fullfile(datadir,files(k).name),'VerTodo');
        V = S.VerTodo;
        assert(all(ismember(V(:,14),1:3)),'%s: cell types must be 1-3',files(k).name);
        n = size(V,1);
        genotype = lower(string(tok{3}));
        parts{k} = table(repmat(genotype + tok{2},n,1),repmat(genotype,n,1),repmat(string(tok{1}),n,1),...
            TYPES(V(:,14))',V(:,8),V(:,7),V(:,13),...
            'VariableNames',{'larva','genotype','side','type','plane','rostroCaudal','lateroMedial'});
    end
    cells = vertcat(parts{:});
    assert(all(ismember(cells.genotype,["wt","mecp2"])),'unexpected genotype labels');

    % rostro-caudal: position along the tectal midline, both hemispheres in the same direction
    right = cells.side == "Rigth";
    cells.rostroCaudal(right) = 1 - cells.rostroCaudal(right);

    % cells beyond the ends of the midline are projected onto its endpoint: no valid position
    ends = cells.rostroCaudal == 0 | cells.rostroCaudal == 1;
    fprintf('Cells projected onto a midline endpoint (excluded from E): %d of %d (%.1f%%)\n',...
        sum(ends),height(cells),100*mean(ends));
    cells.rostroCaudal(ends) = NaN;

    % latero-medial: signed distance to the midline, [-DIST_MAX DIST_MAX] px mapped onto [0 1]
    far = abs(cells.lateroMedial) > DIST_MAX;
    fprintf('Cells farther than %d px from the midline (excluded from F): %d (%.1f%%)\n\n',...
        DIST_MAX,sum(far),100*mean(far));
    cells.lateroMedial = (cells.lateroMedial + DIST_MAX)/(2*DIST_MAX);
    cells.lateroMedial(far) = NaN;

    % larvae imaged bottom-up: flip planes if that matches the other larvae better
    edges = 0:PLANE_MAX+1;
    profile = @(z) histcounts(z,edges,'Normalization','probability')';
    ref = profile(cells.plane(~ismember(cells.larva,FLIPPED)));
    ACTION = ["kept","flipped"];
    disp('Plane orientation check (correlation with the other larvae):')
    for a = FLIPPED(ismember(FLIPPED,cells.larva))
        i = cells.larva == a;
        r = [corr(profile(cells.plane(i)),ref), corr(profile(PLANE_MAX - cells.plane(i)),ref)];
        if r(2) > r(1), cells.plane(i) = PLANE_MAX - cells.plane(i); end
        fprintf('    %s: as stored r=%.2f, flipped r=%.2f -> %s\n',a,r,ACTION((r(2) > r(1)) + 1));
    end
    disp(' ')
end
