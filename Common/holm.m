function q = holm(p)
%%% q = holm(p)
% Holm-Bonferroni adjusted p-values (family-wise error rate); NaN entries are ignored
    q = NaN(size(p));
    ok = find(~isnan(p));
    [ps,order] = sort(p(ok));
    m = numel(ps);
    q(ok(order)) = min(1, cummax(ps(:)' .* (m:-1:1)));
end
