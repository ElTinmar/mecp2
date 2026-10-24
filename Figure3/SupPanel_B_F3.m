function SupPanel_B_F3(...
    N_samp,...
    confusion,...
    decoding,...
    stim_azimuth,...
    ax)

    col = [0 0.447 0.741; 0.85 0.325 0.098];

    nstim = length(stim_azimuth);
    angle_deviation = zeros(nstim,nstim);
    for s1=1:nstim
        for s2 = 1:nstim
            angle_deviation(s1,s2) = abs(stim_azimuth(s1)-stim_azimuth(s2));
        end
    end

    numfish = [size(confusion{1},3) size(confusion{2},3)];
    nrep = size(confusion{1},2);
    ntrials = size(decoding{1}{1},2);
    totstim = nstim * ntrials;

    angular_error = cell(length(N_samp),nrep,2,max(numfish));
    angular_error_wt = [];
    angular_error_mecp2 = [];
    ave_angular_error_wt = [];
    ave_angular_error_mecp2 = [];
    se_angular_error_wt = [];
    se_angular_error_mecp2 = [];
    for s = 1:length(N_samp)
        for rep = 1:nrep
            for fs = 1:2
                for i = 1:numfish(fs)
                    for s1=1:nstim
                        for s2 = 1:nstim
                            angular_error{s,rep,fs,i} = [angular_error{s,rep,fs,i} repmat(angle_deviation(s1,s2),1,confusion{fs}{s,rep,i}(s1,s2))];
                        end
                    end
                end
            end
        end
 
        % Do the average for each stim in totstim over nrep x nfish for mecp2 and wt
        angular_error_wt(s,:) = squeeze(nanmean(reshape([angular_error{s,:,1,:}],totstim,[]),1));  
        angular_error_mecp2(s,:) = squeeze(nanmean(reshape([angular_error{s,:,2,:}],totstim,[]),1));
        % average repetitions within each fish, fish are the unit
        fish_wt{s} = mean(reshape(angular_error_wt(s,:),nrep,[]),1);
        fish_mecp2{s} = mean(reshape(angular_error_mecp2(s,:),nrep,[]),1);
        ave_angular_error_wt(s) = mean(fish_wt{s});
        ave_angular_error_mecp2(s) = mean(fish_mecp2{s});
        se_angular_error_wt(s) = std(fish_wt{s})./sqrt(numfish(1));
        se_angular_error_mecp2(s) = std(fish_mecp2{s})./sqrt(numfish(2));
    end

    axes(ax)
    xx = [N_samp, fliplr(N_samp)];;
    yy = [ave_angular_error_wt+se_angular_error_wt, fliplr(ave_angular_error_wt-se_angular_error_wt)];
    patch(xx,yy,col(1,:),'FaceAlpha',0.5,'EdgeColor','none')
    hold on

    xx = [N_samp, fliplr(N_samp)];;
    yy = [ave_angular_error_mecp2+se_angular_error_mecp2, fliplr(ave_angular_error_mecp2-se_angular_error_mecp2)];
    patch(xx,yy,col(2,:),'FaceAlpha',0.5,'EdgeColor','none')

    plot(N_samp,ave_angular_error_wt,'Color',col(1,:),'LineWidth',2)
    plot(N_samp,ave_angular_error_mecp2,'Color',col(2,:),'LineWidth',2)

    text(0.55, 0.92, 'wild type', 'Units', 'normalized', 'FontSize', 8, 'FontWeight', 'bold', 'Color', col(1,:));
    text(0.55, 0.84, '{\it{}mecp2}-mutant', 'Units', 'normalized', 'FontSize', 8, 'FontWeight', 'bold', 'Color', col(2,:));
    ylabel(['<Absolute angle error [' char(0176) ']>_{fish}'])
    xlabel('Number of neurons')

    p = cellfun(@ranksum, fish_wt, fish_mecp2);
    disp('Angular error, bilateral ranksum test across fish:')
    disp(p)
    disp('Holm-corrected: ')
    disp(holm(p))
end