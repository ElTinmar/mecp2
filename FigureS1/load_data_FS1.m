function cells = load_data_FS1(datafile)
%%% cells = load_data_FS1(datafile)
% one row per segmented cell: larva, genotype, side, type, plane, tectalLine.
% ALL.mat holds a struct array ALL with one element per larva.

    PLANE_MAX = 86;     % last imaging plane
    FLIPPED = [1 2 3];  % wild-type larvae imaged bottom-up

    S = load(datafile,'ALL');
    assert(isfield(S,'ALL'),'no variable ALL in %s',datafile);
    TYPES = ["vglut2a","gad1b","dapi"];  % Color codes 1, 2, 3

    num = @(x) str2double(string(x(:)));  % fields may be stored as numbers or text
    parts = cell(numel(S.ALL),1);
    for k = 1:numel(S.ALL)
        s = S.ALL(k);
        color = num(s.Color);
        n = numel(color);
        assert(all(ismember(color,1:3)),'larva %d: Color codes must be 1-3',k);
        side = string(s.Side(:));
        if isscalar(side), side = repmat(side,n,1); end
        parts{k} = table(...
            repmat(str2double(string(s.Animal(1))),n,1),...
            repmat(lower(string(s.Treatment(1))),n,1),...
            side,...
            TYPES(color)',...
            num(s.Plane),...
            num(s.TectalLine),...
            'VariableNames',{'larva','genotype','side','type','plane','tectalLine'});
    end
    cells = vertcat(parts{:});
    assert(all(ismember(cells.genotype,["wt","mecp2"])),'unexpected genotype labels');

    % both tectal hemispheres on one latero-medial axis (position along the tectal midline)
    right = ismember(cells.side,["Right","Rigth"]);
    cells.tectalLine(right) = 1 - cells.tectalLine(right);

    % cells beyond the ends of the midline are projected onto its endpoint: no valid position
    ends = cells.tectalLine == 0 | cells.tectalLine == 1;
    fprintf('Cells projected onto a midline endpoint (excluded from F): %d of %d (%.1f%%)\n\n',...
        sum(ends),height(cells),100*mean(ends));
    cells.tectalLine(ends) = NaN;

    % larvae imaged bottom-up: flip planes if that matches the other larvae better
    edges = 0:PLANE_MAX+1;
    profile = @(z) histcounts(z,edges,'Normalization','probability')';
    ref = profile(cells.plane(~ismember(cells.larva,FLIPPED)));
    ACTION = ["kept","flipped"];
    disp('Plane orientation check (correlation with the other larvae):')
    for a = FLIPPED
        i = cells.larva == a;
        r = [corr(profile(cells.plane(i)),ref), corr(profile(PLANE_MAX - cells.plane(i)),ref)];
        if r(2) > r(1), cells.plane(i) = PLANE_MAX - cells.plane(i); end
        fprintf('    larva %d: as stored r=%.2f, flipped r=%.2f -> %s\n',a,r,ACTION((r(2) > r(1)) + 1));
    end
    disp(' ')
end
