function panel_C_FS1(counts,ax)

    col = [0 0.447 0.741; 0.85 0.325 0.098];
    wt = counts.genotype == "wt";
    ratio = counts.vglut2a ./ counts.gad1b;  % excitation / inhibition, as in the published figure

    [p,m_wt,se_wt,m_mecp2,se_mecp2] = barplot_ranksum(ax,ratio(wt),ratio(~wt),...
        '{\it{}vglut2a} / {\it{}gad1b}',{'wild type','{\it{}mecp2}-mutant'},col);
    ax.XTickLabelRotation = 0;

    inverse = 1 ./ ratio;
    sem = @(x) std(x)/sqrt(numel(x));
    disp('Excitation / inhibition ratio (vglut2a / gad1b, mean +- SEM)')
    fprintf('WT: %.3f +- %.3f, MeCP2: %.3f +- %.3f\n',m_wt,se_wt,m_mecp2,se_mecp2);
    fprintf('    (gad1b / vglut2a: WT %.3f +- %.3f, MeCP2 %.3f +- %.3f)\n',...
        mean(inverse(wt)),sem(inverse(wt)),mean(inverse(~wt)),sem(inverse(~wt)));
    fprintf('bilateral ranksum test: %.3g\n',p);

    glm = fitglm(table(ratio,categorical(counts.genotype,["wt","mecp2"]),'VariableNames',{'ratio','genotype'}),...
        'ratio ~ genotype','Distribution','gamma');
    fprintf('Gamma GLM (ratio ~ genotype): p=%.3g\n',glm.Coefficients{'genotype_mecp2','pValue'});
    disp(' ')
end
