function define_rois()
files = dir('*_meanFrame.tif');
if ~isempty(files)
    for n = 1:length(files)
        outfile = strrep(files(n).name,'_meanFrame.tif', '_rois.mat');
        if ~exist(outfile,'file')
            img = loadtiff(files(n).name);
            resp = 'n';
            if n>1
                resp = input('use rois from previous file? (y/n)','s');
                if strcmp(resp,'y')
                    prevRoi = load(strrep(files(n-1).name,'_meanFrame.tif', '_rois.mat'));
                    figure(199); clf
                    imshow2(img, [min(img(:)) max(img(:))*0.9], 'InitialMagnification', 'fit')
                    hold on
                    for k = 1:size(prevRoi.rois,2)
                        plot(prevRoi.rois{k}(:,1),prevRoi.rois{k}(:,2))
                    end
                    drawnow
                    resp = input('ok? y to confirm, n to draw new, m to modify','s');
                    val = [0 0];
                    while strcmp(resp,'m')
                        val = str2num(input('value? (x y)','s'));
                        figure(199); clf
                        imshow2(img, [min(img(:)) max(img(:))*0.9], 'InitialMagnification', 'fit')
                        hold on
                        for k = 1:size(prevRoi.rois,2)
                            plot(prevRoi.rois{k}(:,1)+val(1),prevRoi.rois{k}(:,2)+val(2))
                        end
                        drawnow
                        resp = input('ok? y to confirm, n to draw new, m to modify','s');
                    end
                    if strcmp(resp,'y')
                        rois = prevRoi.rois;
                        rois{k}(:,1) = rois{k}(:,1)+val(1);
                        rois{k}(:,2) = rois{k}(:,2)+val(2);
                    end
                end
            end
            if strcmp(resp,'n') || n==1
                rois = getRois(img);
            end
            save(outfile,'rois');
        end
    end
end

    function roi_points = getRois(refimg)
        %based on clicky
        figure(199); clf
        imshow2(refimg, [min(refimg(:)) max(refimg(:))*0.9], 'InitialMagnification', 'fit')
        hold on;
        
        npts = 1;
        colorindex = 0;
        order = get(gca,'ColorOrder');
        nroi = 1;
        while(npts > 0)
            
            figure(199)
            [xv, yv] = (getline(gca, 'closed'));
            if size(xv,1) < 3  % exit loop if only a line is drawn
                break
            end
            
            %draw the bounding polygons and label them
            currcolor = order(1+mod(colorindex,size(order,1)),:);
            plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
            text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
            
            roi_points{nroi} = [xv, yv];
            nroi = nroi + 1;
            colorindex = colorindex+1;
        end
    end
end