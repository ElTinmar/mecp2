function SupMovie_neuralActivity(...
    activity,...
    resp_perStim,...
    ex,...
    stim_subset,...
    stim_duration,...
    thick,...
    speedup,...
    plot_isocurve,...
    videofilename)

    colors = jet(length(stim_subset));

    ind_wt = find(~activity.ismecp2);
    ind_mecp2 = find(activity.ismecp2);

    video_wt = get_video(ind_wt(ex.wt),thick(1));
    video_mecp2 = get_video(ind_mecp2(ex.mecp2),thick(2));

    % account for potential difference in height
    diff_height = size(video_wt,1) - size(video_mecp2,1);
    if diff_height>0
        video_mecp2 = padarray(video_mecp2,[diff_height 0 0 0],0,'post');
    else
        video_wt = padarray(video_wt,[-diff_height 0 0 0],0,'post');
    end

    height = size(video_wt,1);
    cb = cat(3,repmat(linspace(1,0,height)',1,10), zeros(height,10), zeros(height,10)); 
    cb = cat(2, zeros(size(cb)), cb, zeros(size(cb,1),30,3));
    colorscale = repmat(cb,1,1,1,size(video_wt,4));
    final_video = cat(2,video_wt,video_mecp2,colorscale);
   
    export_movie(videofilename,final_video);

    function video = get_video(file,thicken)
        n_cell = size(activity.deltaFoF{file},3);
        n_samples = size(activity.deltaFoF{file},4);
        imrange = prctile(activity.anatomy{file}(:),[1 99]);
        anatomy = repmat(imadjust(activity.anatomy{file},imrange),[1 1 3]);
        height = size(anatomy,1);
        width = size(anatomy,2);
        channels = size(anatomy,3);

        frame_span = round(activity.fps*stim_duration);
        stim_frames = activity.onset:activity.onset+frame_span;

        video = [];

        % one stim at a time
        for idx = 1:length(stim_subset)
            stim = stim_subset(idx);
            norm_dff = prctile(nanmean(activity.deltaFoF{file}(:,stim,:,:)),99,'all');
            for sample = 1:n_samples
                act = zeros(height,width,'single');
                for neuron = 1:n_cell
                    linear_index = activity.cells{file}{neuron};
                    if thicken
                        temp = zeros(height,width);
                        temp(linear_index) = 1;
                        temp = bwmorph(temp,'thicken',1);
                        linear_index = find(temp(:));
                    end
                    act(linear_index) = nanmean(activity.deltaFoF{file}(:,stim,neuron,sample));
                end

                act(act<0) = 0;
                act = act./norm_dff;

                act_rgb = cat(3,...
                    act,...
                    zeros(height,width,'single'),...
                    zeros(height,width,'single'));

                frame = anatomy + act_rgb;
                frame = add_stim(frame,sample,stim_frames,idx);
                if sample == 1
                    [mask_COM, mask_isocurve] = mask_kde(frame,stim,colors(idx,:),file,25);
                end
                if plot_isocurve
                    frame = frame + mask_COM + 0.5*mask_isocurve;
                else
                    frame = frame + mask_COM;
                end
                video = cat(4,video,frame);
            end
        end

        % overlay all together
        frame = anatomy + act_rgb;
        for idx = 1:length(stim_subset)
            [mask_COM, mask_isocurve] = mask_kde(frame,stim_subset(idx),colors(idx,:),file,25);
            frame = frame + mask_COM + 0.5*mask_isocurve;
        end

        for i = 1:20
            video = cat(4,video,frame);
        end

        video(video>1) = 1;
    end

    function [mask_COM, mask_isocurve] = mask_kde(frame,stim,col,file,dotness)
        radius = 3;  
        debug = false;
        
        % use 2D KDE to estimate center of mass of response
        % this is more robust to outliers than mean
        mask_COM = zeros(size(frame));
        mask_isocurve = zeros(size(frame,2),size(frame,1),3);
        [h,w] = meshgrid(1:size(frame,1),1:size(frame,2));
        h_vec = h(:);
        w_vec = w(:); 
        grid = [w_vec h_vec];
  
        respcells = resp_perStim{file,stim};
        positions = zeros(length(respcells),2);
        for idx = 1:length(respcells)
            neuron = respcells(idx);
            linear_index = activity.cells{file}{neuron};
            [y,x] = ind2sub(size(frame),linear_index);
            positions(idx,:) = nanmean([x,y]);
        end
        dens_2D = ksdensity(positions,grid,'BandWidth',[20 20]);
        
        if debug
            figure 
            surf(w,h,reshape(dens_2D,size(frame,2),size(frame,1)))
            drawnow
        end

        md = max(dens_2D);
        [~,loc] = max(dens_2D);
        D = reshape(dens_2D,size(frame,2),[]);
        B = bwboundaries(D >= (md/2));

        for k=1:length(B)
            boundary = B{k};
            
            % create spline
            spl = cscvn(transpose(boundary(:,[1:end 1])));    
            
            % create dotted line
            t = linspace(spl.breaks(1), spl.breaks(end), 5000);
            torm = mod(t,(spl.breaks(end)-spl.breaks(1))/dotness) > ((spl.breaks(end)-spl.breaks(1))/dotness)/2;
            t(torm) = [];

            new_points = fnval(spl, t);
            pixels = round(new_points);

            if debug
                figure
                plot(boundary(:,2),boundary(:,1))
                hold on
                plot(pixels(2,:),pixels(1,:),'o');
                drawnow
            end
                
            for p=1:size(pixels,2)
                mask_isocurve(pixels(1,p),pixels(2,p),:) = col;
            end
        end
        
        centerOfMass = grid(loc,:); 
        mask_COM = insertShape(mask_COM,'Circle',[centerOfMass radius],'Color',col);
        mask_isocurve = permute(mask_isocurve,[2,1,3]);
    end

    function frame = add_stim(frame,sample,stim_frames,stim_num)
        
        outer_radius = 7;
        inner_radius = 6;
        num_stim = length(stim_subset);
        centers = [linspace(10,10+outer_radius*2*num_stim,num_stim);...
            repmat(10,1,num_stim)];
        
        [x,y] = meshgrid(1:size(frame,2),1:size(frame,1));
        x = x(:);
        y = y(:);
        ind_disk = ((x - centers(1,:)).^2 + (y - centers(2,:)).^2) <= outer_radius^2;
        ind_circle = ind_disk ...
            & ((x - centers(1,:)).^2 + (y - centers(2,:)).^2) >= inner_radius^2;
        
        stim_mask_RGB = zeros(size(frame,1),size(frame,2),3);
        for i = 1:num_stim
            stim_mask = zeros(size(frame,1),size(frame,2));
            stim_mask(ind_circle(:,i)) = 1;
            stim_mask_RGB = stim_mask_RGB + pagemtimes(stim_mask,reshape(colors(i,:),[1 1 3]));
        end
        if ismember(sample,stim_frames)
            stim_mask = zeros(size(frame,1),size(frame,2));
            stim_mask(ind_disk(:,stim_num)) = 1;
            stim_mask_RGB = stim_mask_RGB + pagemtimes(stim_mask,reshape(colors(stim_num,:),[1 1 3]));
        end
        frame = frame + stim_mask_RGB;
    end

    function export_movie(filename,video)
        col = [0 0.447 0.741; 0.85 0.325 0.098];
        mov = VideoWriter(filename,'Uncompressed AVI');
        mov.FrameRate = round(speedup*activity.fps);
        mov.open();
        n_frames = size(video,4);
        width = size(video,2); 
        height = size(video,1); 

        legend = zeros(30,200,3);
        legend = insertText(legend, [100 15], '<activity> (AU)',...
            'FontSize',20,...
            'AnchorPoint','Center',...
            'TextColor','w',...
            'BoxOpacity',0);

        for f = 1:n_frames
            img = squeeze(video(:,:,:,f));
            img = insertText(img,[5 10],'wild type',...
                        'FontSize',20,...
                        'TextColor',col(1,:),...
                        'BoxOpacity',0);
            img = insertText(img,[-20+width/2 10],'mecp2-mutant',...
                        'FontSize',20,...
                        'TextColor',col(2,:),...
                        'BoxOpacity',0);
            img = insertText(img,[width 0],'1',...
                    'FontSize',20,...
                    'TextColor','w',...
                    'BoxOpacity',0,...
                    'AnchorPoint','RightTop');
            img = insertText(img,[width height],'0',...
                    'FontSize',20,...
                    'TextColor','w',...
                    'BoxOpacity',0,...
                    'AnchorPoint','RightBottom');
            img(height/2-200/2+1:height/2+200/2,width-30+1:end,:) = flipud(permute(legend,[2,1,3]));
            mov.writeVideo(img);
        end
        mov.close()
    end
end