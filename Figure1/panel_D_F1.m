function panel_D_F1(data,ex,ax)

    cspan = [0 0.2];
    yspan = [-0.05 0.42];
    col = [0 0.447 0.741; 0.85 0.325 0.098];

    % modify the colormap
    b = 7;
    nu = 1.5;
    t = 0.5;
    sigmoid = @(x) 1./(1+exp(-b*(x-t))).^(1/nu);
    n = 256;
    x = 1/n*[1:n];
    cmap = interp1(1:n,parula(n),n*sigmoid(x),'linear','extrap');

    plot_correlations(...
        data.wt.deltaFoF{ex.wt.corr_ind},...
        data.wt.movements{ex.wt.corr_ind},...
        1)
    plot_correlations(...
        data.mecp2.deltaFoF{ex.mecp2.corr_ind},...
        data.mecp2.movements{ex.mecp2.corr_ind},...
        2)

    function plot_correlations(raster,mvt,curr_ax)
        axes(ax(curr_ax));
        cutFrames = find(mvt);
        deltaFoFcut = raster;
        deltaFoFcut(cutFrames,:) = [];
        C = corr(deltaFoFcut);

        % remove the diagonal
        idx_diag = eye(size(C),'logical');
        C(idx_diag) = NaN;

        imagesc(C);
        colormap(ax(curr_ax),cmap)
        caxis(cspan)
        axis image; % that screws colorbar alignment
        xlabel("# Neuron")

        ax(curr_ax).XTick = [];
        ax(curr_ax).YTick = [];

        if curr_ax == 1
            ylabel("# Neuron")
            title('wild type', 'Color',col(1,:))
        else
            title('{\it{}mecp2}-mutant', 'Color',col(2,:))
            c = mycolorbar('Correlation');
            % manually place the colorbar because Matlab is broken with axis image
            c.Position = [0.915 0.347 0.03 0.5325]; 
            c.Label.Position(1) = -0.1;
            c.Label.Position(2) = 0.105;
            c.Label.Color = 'k';
            c.Label.FontSize = 8;
        end

        
    end
end