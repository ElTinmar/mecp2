function panel_B_FS1(counts,ax)

    col = [0 0.447 0.741; 0.85 0.325 0.098];
    wt = counts.genotype == "wt";
    TYPES = ["vglut2a","gad1b"];

    % rank-sum tests between genotypes, Holm-corrected over the two cell types
    disp('Number of cells per larva (x10^3, mean +- SEM)')
    Y = cell(2,2);  % cell type x genotype
    p_geno = NaN(1,2);
    for k = 1:2
        y = counts.(char(TYPES(k)))/1000;
        Y(k,:) = {y(wt), y(~wt)};
        p_geno(k) = ranksum(y(wt),y(~wt));
        fprintf('%s  WT: %.2f +- %.2f, MeCP2: %.2f +- %.2f\n',TYPES(k),...
            mean(y(wt)),sem(y(wt)),mean(y(~wt)),sem(y(~wt)));
    end
    p_geno = holm(p_geno);
    fprintf('genotype, bilateral ranksum test, Holm-corrected: vglut2a %.3g, gad1b %.3g\n\n',p_geno);

    % Poisson GLME on counts, larva as random effect, overdispersion estimated
    n = height(counts);
    T = table(categorical([counts.larva; counts.larva]),...
        categorical([counts.genotype; counts.genotype],["wt","mecp2"]),...
        categorical([repmat("vglut2a",n,1); repmat("gad1b",n,1)],TYPES),...
        [counts.vglut2a; counts.gad1b],...
        'VariableNames',{'larva','genotype','type','count'});
    glme = fitglme(T,'count ~ type*genotype + (1|larva)',...
        'Distribution','Poisson','Link','log','DispersionFlag',true);
    disp('GLME: count ~ type*genotype + (1|larva), Poisson, log link, dispersion estimated')
    disp(glme.Coefficients)
    names = string(glme.CoefficientNames);
    g = names == "genotype_mecp2";
    t = names == "type_gad1b";
    gt = contains(names,":");
    contrasts = {'genotype, vglut2a', g; 'genotype, gad1b', g | gt;...
                 'gad1b vs vglut2a, wild type', t; 'gad1b vs vglut2a, mecp2-mutant', t | gt};
    p_glme = NaN(1,4);
    for k = 1:4
        [p_glme(k),F,df1,df2] = coefTest(glme,double(contrasts{k,2}));
        fprintf('    %s: F(%g,%g)=%.2f, p=%.3g\n',contrasts{k,1},df1,df2,F,p_glme(k));
    end
    disp(' ')

    % plot: bars = mean, error bars = SEM, dots = larvae
    x = [1 2; 3.5 4.5];  % cell type x genotype
    axes(ax)
    hold on
    for k = 1:2
        for gi = 1:2
            y = Y{k,gi};
            m = mean(y);
            bar(x(k,gi),m,0.9,'FaceColor',col(gi,:),'EdgeColor','k','LineWidth',1);
            plot([x(k,gi) x(k,gi)],m + [0 sem(y)],'k','LineWidth',1);
            plot(x(k,gi) + [-.1 .1],[1 1]*(m + sem(y)),'k','LineWidth',1);
            swarmchart(repmat(x(k,gi),numel(y),1),y,'.k');
        end
    end
    top = max(cellfun(@max,Y),[],'all');
    xlim([0.3 5.2])
    ylim([0 1.6*top])
    bracket(x(1,:),1.06*top,p_geno(1));
    bracket(x(2,:),1.06*top,p_geno(2));
    bracket(x(:,1)',1.24*top,p_glme(3));
    bracket(x(:,2)',1.42*top,p_glme(4));
    ax.XTick = mean(x,2);
    ax.XTickLabel = {'{\it{}vglut2a}','{\it{}gad1b}'};
    ylabel('Cell number (\times10^3)')

    function bracket(xx,y,p)
        plot(xx([1 1 2 2]),y - [0.02 0 0 0.02]*diff(ylim),'k','LineWidth',0.75);
        text(mean(xx),y,stars(p),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
    end
end

function s = sem(x)
    s = std(x)/sqrt(numel(x));
end

function s = stars(p)
    s = 'ns';
    if p < 0.001, s = '***'; elseif p < 0.01, s = '**'; elseif p < 0.05, s = '*'; end
end
