% found using plots of ksdensity
edges = 0:0.01:6;
d = diff(edges)/2;
centers = edges(1:end-1)+d;
col = [0 0.447 0.741; 0.85 0.325 0.098];

snr_wt = get_snr(neuralActivity.wt,edges);
snr_mecp2 = get_snr(neuralActivity.mecp2,edges);

figure
p = semilogy(centers,mean(snr_wt));
p.Marker = '.';
p.Color = col(1,:);
p.LineStyle = 'none';
p.MarkerSize = 5;
hold on

p1 = semilogy(centers,mean(snr_mecp2));
p1.Marker = '.';
p1.Color = col(2,:);
p1.LineStyle = 'none';
p1.MarkerSize = 5;

function dist_snr = get_snr(data,edges)

    for file = 1:numel(data.deltaFoF_ori)
        DFoF = data.deltaFoF_ori{file};
        n_cells = size(DFoF,2);
        n_samples = size(DFoF,1);
        sigma = zeros(n_cells,1);
        for c = 1:n_cells
            [d,x] = ksdensity(DFoF(:,c));
            [~,indPeak] = max(d);
            xFit=x(1:indPeak);
            dataToFit=d(1:indPeak)/n_samples;
            [sigma(c),~,~] = mygaussfit(xFit',dataToFit);
        end

        snr = prctile(DFoF,95)./sigma';
        dist_snr(file,:) = histcounts(snr,edges,'Normalization','probability');
    end
end