function counts = count_cells_FS1(cells)
%%% counts = count_cells_FS1(cells)
% number of cells of each type per larva

    [larva,first,idx] = unique(cells.larva);
    n = @(type) accumarray(idx,double(cells.type == type));
    counts = table(larva,cells.genotype(first),n("vglut2a"),n("gad1b"),n("dapi"),...
        'VariableNames',{'larva','genotype','vglut2a','gad1b','dapi'});

    fprintf('WT: %d larvae, MeCP2: %d larvae\n',sum(counts.genotype=="wt"),sum(counts.genotype=="mecp2"));
    disp(counts)
end
