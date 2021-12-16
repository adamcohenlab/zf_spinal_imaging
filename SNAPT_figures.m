cells = [998 1869 2200 2310 2453 2457 2463];
% cells = [2310 2453];

fits = {};
iter = 1;
for cs = cells
    mov = loadtiff(['..\Data\cell_' num2str(cs) '.tif']);
    img = loadtiff(['..\Data\corr_cell_' num2str(cs) '.tif']);
    
    %% interpolate spike mov to correct for rolling shutter
    [ySize, xSize, ~] = size(mov);
    srImg = 996.3;
    xAxImg = 1/srImg:1/srImg:size(mov,3)/srImg;
    center = floor(ySize/2);
    % dtimgCorr([center, center],:) = dtimg([center, center],:);
    spikemovCorr = zeros(size(mov));
    for n = 1:center
        pixelMat = mov([center + n-1, center - (n-1)],:,:);
        pixelMat = reshape(pixelMat, 2*size(pixelMat,2), size(pixelMat,3))';
        xq = xAxImg - 9.74436e-6 * (n-1);
        pixelInterp = interp1(xAxImg, pixelMat, xq, 'cubic');
        pixelInterp = reshape(pixelInterp', 2, size(mov, 2), size(mov, 3));
        spikemovCorr([center + n-1, center - (n-1)],:,:) = pixelInterp;
    end
    spikemov = spikemovCorr(:,:,2:end);
    
    % for n = 1:size(spikemov,3)
    %     spikemov(:,:,n) = imgaussfilt(spikemov(:,:,n),3);
    % end
    %%
    background = movmean(quantile(mean(spikemov(:,:,1:10),3), .2),80);
    movbg = spikemov - 100;
    % movbg = spikemov - background;
    dfmov = movbg-mean(movbg(:,:,1:3),3);
    dffmov = dfmov./mean(movbg(:,:,1:3),3);
    %%
    files = dir(['..\Data\' num2str(cs) '_Path*.roi']);
    
    for rs = 1:length(files)
        tmpr = ReadImageJROI([files(rs).folder '\' files(rs).name]);
        X = tmpr.mnCoordinates(:,1);
        Y = tmpr.mnCoordinates(:,2);
        %% move roi along line and calculate average, based on polyLineKymo and apply_clicky
        dX = diff(X); % Vector pointing from start to end of line
        dY = diff(Y);
        nSeg = length(dX); % number of line segments
        
        dX0 = dX./(dX.^2 + dY.^2).^.5; % Unit vector along line
        dY0 = dY./(dX.^2 + dY.^2).^.5;
        
        dXp = -dY0; % Unit vector perpendicular to line
        dYp = dX0;
        
        L = (dX.^2 + dY.^2).^.5;  % Length of the line
        c = 1;
        dL = 1;
        dP = 8;
        tmproi = {};
        for k = 1:nSeg
            lSteps = 0:dL:L(k);  % displacements along the line
            lSteps = lSteps(1:end-1);
            nSteps = length(lSteps);  % total number of steps
            
            % Define the coordinates of all the rectangles
            for j = 1:nSteps
                Xc = X(k) + lSteps(j)*dX0(k); % box center
                Yc = Y(k) + lSteps(j)*dY0(k);
                Xs = [Xc - dX0(k)*dL/2 - dXp(k)*dP/2; % Make closed paths (5 points)
                    Xc - dX0(k)*dL/2 + dXp(k)*dP/2;
                    Xc + dX0(k)*dL/2 + dXp(k)*dP/2;
                    Xc + dX0(k)*dL/2 - dXp(k)*dP/2;
                    Xc - dX0(k)*dL/2 - dXp(k)*dP/2];
                Ys = [Yc - dY0(k)*dL/2 - dYp(k)*dP/2;
                    Yc - dY0(k)*dL/2 + dYp(k)*dP/2;
                    Yc + dY0(k)*dL/2 + dYp(k)*dP/2;
                    Yc + dY0(k)*dL/2 - dYp(k)*dP/2;
                    Yc - dY0(k)*dL/2 - dYp(k)*dP/2];
                Xs(Xs > xSize) = xSize; Xs(Xs < 1) = 1;
                Ys(Ys > ySize) = ySize; Ys(Ys < 1) = 1;
                tmproi{c} = [Xs, Ys];
                c = c + 1;
            end
        end
        [~,nroi] = size(tmproi);
        dffvec = tovec(dffmov);
        imgvec = tovec(img);
        [x, y] = meshgrid(1:xSize, 1:ySize);
        itrace = zeros(size(dffmov,3), nroi);
        lineintens = zeros(nroi,1);
        for j = 1:nroi
            xv = tmproi{j}(:,1);
            yv = tmproi{j}(:,2);
            inpoly = inpolygon(x,y,xv,yv);
            itrace(:,j) = mean(dffvec(inpoly,:))';
            lineintens(j) = mean(imgvec(inpoly,:));
        end
        %%
        % smooth in space
        itraceav = movmean(itrace,40,2);
        lineintens = movmean(lineintens,40);
        %smooth in time
        itraceav = movmean(itraceav,1,1);
        % scale px/um
        pxscale = 2.13;
        pxAx = (1:nSeg)/pxscale;
        figure
        nexttile
        plot(pxAx, lineintens)
        ylim([-.2 .8])
        yline(.1)
        nexttile
        
        % kymograph
        
        
        imagesc(10:30, pxAx, nrm(itraceav(10:30,:))')
        cmocean('thermal')
        xlabel('Time (ms)')
        ylabel('Distance (\mum)')
        hold on
        
        % get frame of max value at each positio in space
        [~, ma] = max(itraceav);
        scatter(ma, pxAx, '.k')
        
        % framerate in Hz
        srImg = 996.3;
        xv = pxAx;
        yv = ma;
%         % fit line to peak position
%         fobj = fit(xv(:), yv(:), 'poly1');
%         slope = fobj.p1;
%         intercept = fobj.p2;
%         
%         % conduction velocity peak
%         cv = ((1/slope)*srImg)/(1e6);
        % disp(['cv = ' num2str(cv) ' m/s'])
        % plot([fobj(1) fobj(size(itraceav,2))], [1 nSeg], 'k')
        halfmaxtime = zeros(nSeg,1);
        maxtime = zeros(nSeg,1);
        downhalfmaxtime = zeros(nSeg,1);
        % fit rising edge
        smoothedFits = zeros(length(1:.1:40),nSeg);
        for n = 1:nSeg
            [~, mi] = min(itraceav(1:ma(n), n),[],1);
            % calculate threshol
            thr = itraceav(mi,n) + 0.3*(itraceav(ma(n),n)-itraceav(mi,n));
            % find last frame that crosses threshold before max
            mi2 = find(diff(itraceav(mi:ma(n),n)>thr)==1,1,'last');
            %     mi2 = find(itraceav(mi:ma(n),n)>thr,1)-1;
            xf = (mi+mi2):ma(n);
            yf = itraceav((mi+mi2):ma(n),n);
            % fit smoothing spline
            opts = fitoptions( 'Method', 'SmoothingSpline' );
            opts.SmoothingParam = 0.6;
            xv = 1:size(itraceav,1);
            yv = itraceav(:,n);
            % fit spline
            fobj = fit(xv(:), yv(:), 'smoothingspline', opts);
            % evaluate fit and keep interpolated traces for later
            smoothedFits(:,n) = fobj(1:.1:40);
            % find max by finding zeros of 1st derrivative
            tmpextr = fnzeros(fnder(fobj.p));
            % keep value of maximum that is negative in 2nd derrivative
            tmpidma = fnval(fnder(fnder(fobj.p)),tmpextr)<0;
            % keep value of minimum that is positive in 2nd derrivative
            tmpidmi = fnval(fnder(fnder(fobj.p)),tmpextr)>0;
            tmpmax = tmpextr(1,tmpidma(1,:));
            tmpmin = tmpextr(1,tmpidmi(1,:));
            % evaluate function at max positions
            tmpmaxval = fobj(tmpmax);
            % get index of highest value
            [~, tmpmaxid] = max(tmpmaxval);
            maxtime(n) = tmpmax(tmpmaxid);
            % keep only minima after the maximum
            tmpmin = tmpmin(tmpmin>maxtime(n));
            % if there is no minimum use the last frame
            if ~isempty(tmpmin)
                % evaluate function at min position after peak
                tmpminval = fobj(tmpmin);
                % get index of lowest value after the peak
                [~, tmpminid] = min(tmpminval);
                % get time of lowest minimum
                mintime = tmpmin(tmpminid);
            else
                mintime = xv(end);
            end
            % get the time of the half max after the peak
            downhalfmaxval = 0.5*(fobj(maxtime(n)) - fobj(mintime)) + fobj(mintime);
%             tmptime = slmsolve(fobj.p, downhalfmaxval);
%             tmpid = find(tmptime>maxtime(n),1);
%             downhalfmaxtime(n) = tmptime(tmpid);
            downhalfmaxtime(n) = 0;
            
            if length(xf)>1
                % perform fit
                fobj = fit(xf(:), yf(:), 'poly1');
                % slope
                a = fobj.p1;
                % intercept
                b = fobj.p2;
                %calculate threshold crossing
                thrcross = thr/a - b/a;
                % calculate interpolated max time
%                 maxtime(n) = itraceav(ma(n), n)/a - b/a;
                % calculate halfway between threshold crossing and max
                halfmaxtime(n) = thrcross + (maxtime(n)-thrcross)/2;
            else
                halfmaxtime(n) = NaN;
            end
        end
        scatter(downhalfmaxtime(lineintens>.1), pxAx(lineintens>.1), '.w')
        fits(iter).smoothed = smoothedFits;
        fits(iter).lineintens = lineintens;
        fits(iter).cellnum = cs;
        fits(iter).halfmax = halfmaxtime;
        fits(iter).downhalfmax = downhalfmaxtime;
        fits(iter).pxAx = pxAx;
        % fit halfmax time
%         fobj = fit(pxAx(:), halfmaxtime(:), 'poly1');
%         fits(iter).halfmaxfit = fobj;
%         slope = fobj.p1;
%         intercept = fobj.p2;
        % conduction velocity halfmax
%         cv1 = ((1/slope)*srImg)/(1e6);
%         fits(iter).cv1 = cv1;
%         plot(fobj(pxAx), pxAx, 'w')
        
        scatter(maxtime(lineintens>.1), pxAx(lineintens>.1), '.g')
        fits(iter).maxtime = maxtime;
        % fit max time
%         fobj = fit(pxAx(:), maxtime(:), 'poly1');
        fobj = fit(pxAx(lineintens>.1)', maxtime(lineintens>.1), 'poly1');
        fits(iter).maxtfit = fobj;
        slope = fobj.p1;
        intercept = fobj.p2;
        % conduction velocity max
        cv2 = ((1/slope)*srImg)/(1e6);
        fits(iter).cv2 = cv2;
        plot(fobj(pxAx), pxAx, 'g')
        
        title(['cv2 = ' num2str(cv2) ' m/s'])
        
        cm = cmocean('phase', 400);
        nexttile
        stairs(.5:1:39.5,nrm(itraceav(:,[1 round(size(itrace,2)/2) size(itrace,2)])))
        set(gca, 'ColorOrder', cm([1 160 320],:))
        xlabel('Time (ms)')
        ylabel('Normalized fluorescence')
        nexttile
        hold on
        %plot a few smoothed traces
        plot(1:.1:40, nrm(fits(iter).smoothed(:,[1 round(size(itrace,2)/2) size(itrace,2)])))
        set(gca, 'ColorOrder', cm([1 160 320],:))
        nexttile
        plot(fits(iter).pxAx(lineintens>.1), (fits(iter).downhalfmax(lineintens>.1) - fits(iter).maxtime(lineintens>.1)),'k')
        
%         saveas(gcf, ['X:\Lab\Labmembers\Urs Boehm\Data\7-zArchon\20210812_STA_movies\' num2str(cs) '_' num2str(rs) '.fig'])
        iter = iter + 1;
    end
end

%%
for cs = 1:7
    %%
    img = loadtiff(['..\Data\corr_cell_' num2str(cells(cs)) '.tif']);
    files = dir(['..\Data\' num2str(cells(cs)) '_Path*.roi']);
    dtimg = nan(ySize, xSize);
    fitids = find([fits.cellnum] == cells(cs));
    iter = 1;
    for rs = 1:length(files)
        tmpr = ReadImageJROI([files(rs).folder '\' files(rs).name]);
        X = tmpr.mnCoordinates(:,1);
        Y = tmpr.mnCoordinates(:,2);
        %%
        %     mask = zeros(ySize, xSize);
        %     mask(sub2ind([ySize xSize],Y,X)) = 1;
        %     figure
        %     imagesc(imgaussfilt(mask,4).*img)
        
        
        [x, y] = meshgrid(1:xSize, 1:ySize);
        
        for n = 1:length(fits(fitids(rs)).pxAx)
            if fits(fitids(rs)).lineintens(n) >.1
                tmpimg = nan(ySize, xSize);
                xv = sin(0:.1:2*pi).*20 + X(n+1);
                yv = cos(0:.1:2*pi).*20 + Y(n+1);
                inpoly = inpolygon(x, y, xv, yv);
                tmpimg(inpoly) = fits(fitids(rs)).maxtime(n)-fits(fitids(rs)).maxtime(1);
                dtimg(:,:,iter) = tmpimg;
                iter = iter + 1;
            end
        end
    end
    dtimg = mean(dtimg,3,'omit');
    figure
    imagesc(dtimg)
    %%
    subframeT = 0.05; % ms
    initialT = -1; % ms
    finalT = 3; % ms
    sigma = .4; % ms.  This is how long the flash lasts at each pixel.  Set approximately ~err/2.  Adjust to get a reasonable-looking movie.
    
    times = initialT:subframeT:finalT;
    nSNAPT = length(times);
    
    GaussPeaksmov = zeros(ySize,xSize,nSNAPT);
    for q = 1:nSNAPT
        %  GaussPeaksmov(:,:,q) = exp(-(dtimg-times(q)*ones(ysize,xsize)).^2/(2*sigma^2));
        GaussPeaksmov(:,:,q) = exp(-(times(q)-dtimg).^2./(2.*(sigma).^2));
    end
    %%
    ampimg = imgaussfilt(img,3);
    ampimg = mat2gray(ampimg, double([0 max(ampimg(:))]));
    superlocmov = GaussPeaksmov.*repmat(ampimg, [1 1 nSNAPT]);  % Scale pixels by the appropriate amplitude
    % superlocmov = GaussPeaksmov.*repmat(goodpix, [1 1 nSNAPT]);  % Scale pixels by the appropriate amplitude
    superlocmov(isnan(superlocmov)) = 0;
    moviefixsc(superlocmov)
    %%
    % Make a pretty color movie
    superlocColormov = zeros(ySize, xSize, 3, nSNAPT);
    colorimg3 = repmat(ampimg, [1 1 3]);
    
    cmin4 = 0;
    cmax4 = 0.4;
    for j = 1:length(times)
        superlocColormov(:,:,:,j) = .7*colorimg3 + .3*grs2rgb(superlocmov(:,:,j), colormap(hot), cmin4, cmax4);
    end
    
    figure(20)
    % while(1)
    clear M
    for j = 1:nSNAPT
        imshow2(superlocColormov(:,:,:,j), 'InitialMagnification', 100)
        text(135,15,[num2str(times(j), '%+3.2f') ' ms'], 'FontSize', 10, 'color', [1 1 1], 'HorizontalAlignment', 'right')
        drawnow
        M(j) = getframe(gca);
        %     pause(0.1)
        % end;
    end
%     v = VideoWriter(['X:\Lab\Labmembers\Urs Boehm\Data\7-zArchon\20210812_STA_movies\neurite_mov_spline' num2str(cells(cs)) '.avi'], 'Uncompressed AVI');
%     open(v)
%     writeVideo(v,M)
%     close(v)
    
%     % make time coded static image
%  nanvals = isnan(dtimg);
% cm = cmocean('phase', 400);
% colimg = grs2rgb(dtimg, cm(1:.8*400,:), -1, 2);
% for n = 1:3
%     tmpimg = colimg(:,:,n);
%     tmpimg(nanvals) = 1;
%     colimg(:,:,n) = tmpimg;
% end
% hsvimg = rgb2hsv(colimg);
% hsvimg(:,:,3) = hsvimg(:,:,3).*ampimg;
% hsvimg1 = hsvimg;
% hsvimg(:,:,1) = mod(hsvimg(:,:,1)+.5,1);
% hsvimg2 = hsvimg;
% 
% figure
% nexttile
% imshow(hsv2rgb(hsvimg1))
% nexttile
% imshow(imcomplement(hsv2rgb(hsvimg2)))
% nexttile
% imagesc([min(dtimg(:)) max(dtimg(:))], [0 1], linspace(min(dtimg(:)),max(dtimg(:))))
% colormap(cm(1:.8*400,:))
%     saveas(gcf, ['X:\Lab\Labmembers\Urs Boehm\Data\7-zArchon\20210812_STA_movies\' num2str(cells(cs)) 'time_coded.fig'])
%     saveastiff(uint8(imcomplement(hsv2rgb(hsvimg2)).*255), ['X:\Lab\Labmembers\Urs Boehm\Data\7-zArchon\20210812_STA_movies\' num2str(cells(cs)) '_time_white.tif'])
end
%% make time coded static image
nanvals = isnan(dtimg);
cm = cmocean('phase', 400);
colimg = grs2rgb(dtimg, cm(1:.8*400,:));
codedimg = repmat(ampimg,1, 1, 3);
for n = 1:3
    tmpimg = colimg(:,:,n);
    tmpimg(nanvals) = 1;
    colimg(:,:,n) = tmpimg;
    codedimg(:,:,n) = codedimg(:,:,n).*colimg(:,:,n);
end
figure
nexttile
imshow(colimg)
nexttile
imshow(codedimg)
%% make time coded static image in HSV space
nanvals = isnan(dtimg);
cm = cmocean('phase', 400);
colimg = grs2rgb(dtimg, cm(1:.8*400,:));
for n = 1:3
    tmpimg = colimg(:,:,n);
    tmpimg(nanvals) = 1;
    colimg(:,:,n) = tmpimg;
end
hsvimg = rgb2hsv(colimg);
hsvimg(:,:,3) = hsvimg(:,:,3).*ampimg;
hsvimg(:,:,1) = mod(hsvimg(:,:,1)+.5,1);

codedimg = hsv2rgb(hsvimg);

figure
nexttile
imshow(colimg)
nexttile
imshow(codedimg)
nexttile
imshow(imcomplement(codedimg))