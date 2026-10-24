function [p,tbl] = mixed_anova(Y1,Y2)
%%% [p,tbl] = mixed_anova(Y1,Y2)
% mixed-design ANOVA: genotype between fish (Y1, Y2: fish x levels), level within fish.
% Fish with a missing level are excluded. p = [genotype, level, interaction],
% Greenhouse-Geisser corrected; tbl adds partial (pEta2) and generalized (gEta2) eta squared.
% tbl.Properties.Description lists the fish included and excluded, and the Greenhouse-Geisser epsilon.
    Y = [Y1; Y2];
    genotype = categorical([zeros(size(Y1,1),1); ones(size(Y2,1),1)]);
    keep = ~any(isnan(Y),2);
    n1 = size(Y1,1);
    out = find(~keep)';
    summary = sprintf('Fish in ANOVA: %d WT, %d MeCP2', sum(keep(1:n1)), sum(keep(n1+1:end)));
    if ~isempty(out)
        summary = sprintf('%s; excluded (missing levels): WT #[%s], MeCP2 #[%s]', summary, ...
            num2str(out(out<=n1)), num2str(out(out>n1)-n1));
    end
    t = array2table(Y(keep,:));
    t.genotype = genotype(keep);
    within = table(categorical((1:size(Y,2))'),'VariableNames',{'stim'});
    rm = fitrm(t,sprintf('Var1-Var%d ~ genotype',size(Y,2)),'WithinDesign',within);
    tbl = ranova(rm,'WithinModel','stim');
    tbl.Properties.Description = sprintf('%s; Greenhouse-Geisser epsilon: %.4f', summary, epsilon(rm).GreenhouseGeisser);
    terms = {'genotype','(Intercept):stim','genotype:stim'};
    p = tbl{terms,'pValueGG'}';
    ss = tbl{terms,'SumSq'};
    ss_err = [tbl{'Error','SumSq'}; tbl{'Error(stim)','SumSq'}; tbl{'Error(stim)','SumSq'}];
    tbl.pEta2 = NaN(height(tbl),1);
    tbl.gEta2 = NaN(height(tbl),1);
    tbl{terms,'pEta2'} = ss./(ss+ss_err);
    tbl{terms,'gEta2'} = ss./(ss+tbl{'Error','SumSq'}+tbl{'Error(stim)','SumSq'});
end
