function panel_A_F1(data,ex,ax)
%% TODO fix the size of the ticks 

    cspan = [0 0.75];
    yspan = [-0.05 0.42];
    col = [0 0.447 0.741; 0.85 0.325 0.098];
    fps = 1000/256;

    % high pass filter (1/cutoff_frequency ~= 30s)
    L = 101;
    apf = zeros(1,L);
    apf(ceil(L/2)) = 1;
    lpf = 1/L*ones(1,L);
    hpf = apf - lpf;

    plot_raster('wt')
    plot_raster('mecp2')

    function plot_raster(cond)

        % get the data
        R = data.(cond).raster{ex.(cond).raster_ind};
        count = sum(R,2)./size(R,2);
        frames = ex.(cond).raster_frames;

        X = data.(cond).deltaFoF_ori{ex.(cond).raster_ind}';
        X = X(:,frames);

        % apply high pass filter 
        X = conv2(X,hpf,'same');
        X = X - prctile(X(:),80);
        
        % plot raster
        if strcmp(cond,'wt')
            add = 0;
        else
            add = 2;
        end

        axes(ax(1+add));
        imagesc(X);
        ax(1+add).XTick = 1:3*60*fps:length(frames);
        ax(1+add).XTickLabel = [];
        caxis(cspan)
        brighten(-1)
        colormap(1-gray)
        c = mycolorbar('\Delta{}F/F');
        c.Position(1) = 0.96;
        c.Label.Position(1) = -0.09;
        c.Label.Position(2) = 0.375;
        c.Label.Color = 'w';
        c.Label.FontSize = 8;
        if strcmp(cond,'wt')
            ylabel({'{\bf{}wild type}','\color{black}# Neuron'},'Color',col(1,:));
        else
            ylabel({'{\bf{\it{}mecp2}-mutant}','\color{black}# Neuron'}, 'Color',col(2,:));
        end
        
        % plot fraction of active cells
        axes(ax(2+add));
        b=bar(count(frames)*100);
        b.BarWidth=1;
        
        xlim([1 length(frames)])
        ylim([0 10])
        ax(2+add).XTick = 1:3*60*fps:length(frames);
        ax(2+add).XTickLabel = [];
        ax(2+add).YTick = [0:4:8];

        % align with raster
        ax(2+add).Position(3) = ax(1+add).Position(3);

        ylabel({'% Active','neuron'})
        if strcmp(cond,'mecp2')
            xlabel('Time (min)')
            ax(2+add).XTickLabel = 0:3:15;
            b.FaceColor=col(2,:);
            b.EdgeColor=col(2,:);
        else
            b.FaceColor=col(1,:);
            b.EdgeColor=col(1,:);
        end
    end
end