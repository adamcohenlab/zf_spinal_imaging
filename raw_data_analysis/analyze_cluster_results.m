folder = dir();
for n = 1:length(folder)
    try
    if folder(n).isdir
        file = dir([folder(n).name '\*.Vp']);
        if ~isempty(file)
            % open vnr data
            [dat,srate] = import2pdaq([folder(n).name '\'], file.name,'a');
            % bandpass between 300 and 100 Hz to filter our signal
            V = bandpass(dat-mean(dat),[300 1000],srate);
            %hard coded image parameters
            srImg = 996.3;
            %image time axis
            xAxImg = 1/srImg:1/srImg:25000/srImg;
            % vnr time axis
            xAxVNR = 1/srate:1/srate:length(V)/srate;
            
            
            stdV = movstd(V,srate*0.01);        % 10ms sliding window std
            stdV = movmean(stdV,srate*0.01);
            [N,edges] = histcounts(stdV,200);   % histogram of std trace
            %     [~, maxIdx] = max(N);                       % histogram peak
            %     noise = [N(1:maxIdx) fliplr(N(1:maxIdx-1))];% mirror left side of histogram
            %     fobj = fit((1:length(noise))', noise', 'gauss1');   % fit gausian
            fobj = fit((1:length(N))', N', 'gauss1');   % fit gausian
            thr = edges(round(3*fobj.c1/sqrt(2) +fobj.b1));                          % peak threshold is 4*std of gaussian fit
            [peaks,peakInd] = findpeaks(stdV,'MinPeakHeight',thr, 'MinPeakWidth', srate*0.01);
            bursts = peakInd;
            % instantaneous frequency of each burst
            instFreq = [0; 1./diff(peakInd./srate)];
            % get rid of empty data
            if isempty(instFreq)
                instFreq = [];
            end
            if isempty(instFreq)
                instFreq = [];
            end
            boutIdx = 0;
            boutStart = [];
            boutStartIdx = [];
            % vector true when bursts occur
            tmp = stdV-thr;
            tmp(tmp<0) = 0;
            for k = 1:length(bursts)
                % get burst start
                burstStart = find(tmp(1:peakInd(k))==0,1, 'last');
                burstEnd = find(tmp(peakInd(k):end)==0,1, 'first');
                burstEnd = burstEnd + peakInd(k);
                if instFreq(k)<5
                    % get bour start
                    boutIdx = boutIdx + 1;
                    boutStart(boutIdx) = bursts(k);
                    boutStartIdx(boutIdx) = k;
                end
            end
            boutEndIdx = [];
            boutEnd = [];
            for k = 1:length(boutStartIdx)
                if k < length(boutStartIdx)
                    % get burst end
                    boutEndIdx(k) = boutStartIdx(k+1) - 1;
                    boutEnd(k) = bursts(boutEndIdx(k));
                else
                    boutEndIdx(k) = length(bursts);
                    boutEnd(k) = bursts(boutEndIdx(k));
                end
            end
            
            % clean up bursts that are not part of a swim bout
            noBurst = (boutEnd - boutStart) == 0;
            burstsStartTmp = zeros(size(bursts));
            burstsStartTmp(boutStartIdx) = 1;
            burstsEndTmp = zeros(size(bursts));
            burstsEndTmp(boutEndIdx) = 1;
            bursts(boutStartIdx(noBurst)) = [];
            instFreq(boutStartIdx(noBurst)) = [];
            burstsStartTmp(boutStartIdx(noBurst)) = [];
            boutStartIdx = find(burstsStartTmp);
            boutStart = bursts(boutStartIdx);
            burstsEndTmp(boutEndIdx(noBurst)) = [];
            boutEndIdx = find(burstsEndTmp);
            boutEnd = bursts(boutEndIdx);
            
            figure
            f1 = subplot(3,1,1);
            plot(xAxVNR,V)
            ylabel('V/mV')
            f2 = subplot(3,1,2);
            plot(xAxVNR,stdV)
            hold on
            plot(xAxVNR(bursts),stdV(bursts),'.k');
            plot(xAxVNR(boutStart), stdV(boutStart),'g.')
            plot(xAxVNR(boutEnd), stdV(boutEnd),'r.')
            ylabel('std 10ms window /a.u.')
            f3 = subplot(3,1,3);
            plot(xAxVNR(bursts),instFreq,'o');
            ylabel('inst. frequency / Hz')
            xlabel('time/s') 
            linkaxes([f1 f2 f3], 'x')
            
            clear burstImg
            burstTime = xAxVNR(bursts);
            for k = 1:length(burstTime)
                [~, burstImg(k)] = min((xAxImg-burstTime(k)).^2);
            end
            save(strrep(file.name,'.Vp', '_VNR.mat'),'V', 'dat', 'bursts', 'burstImg','stdV','boutStart',...
                'boutStartIdx', 'boutEnd', 'boutEndIdx', 'instFreq');
        end
    end
    catch
        disp(['error in file ' folder(n).name '\' file.name])
    end
end
%%
files = dir('*meanFrame.tif');


for k = 1:length(files)
    fileN = str2num(cell2mat(regexp(files(k).name,'\d*','match')));
    if ~isempty(dir([ '*tr' num2str(fileN) '*_VNR.mat']))
        img = loadtiff(files(k).name);
        vnrFile = dir([ '*tr' num2str(fileN) '*_VNR.mat']);
        load(vnrFile.name);
        load(strrep(files(k).name,'meanFrame.tif','ica_pca_selected.mat'));
        load(strrep(files(k).name,'meanFrame.tif','rois.mat'));
%         load(strrep(files(k).name,'meanFrame.tif','ica_pca_trace.mat'));
%         tr = [];
%         for n = 1:length(trace)
%             maxFrame = length(trace{n}(:,1));
%             if ~isempty(kept{n})
%                 tmp = trace{n}(:,kept{n});
%                 
%                 tmp = correct_bleaching(tmp(200:maxFrame-round(maxFrame/100)));
%                 tmp = tmp(301:end);
%                 tmp = (tmp-mean(tmp))./mean(tmp);
%             else
%                 tmp = zeros(size(trace{n}(500:maxFrame-round(maxFrame/100),1)));
%             end
%             tr(:,n) = tmp;
%         end
        a = rois;
        instFreq = 1./diff(bursts./srate);
        timf = zeros(size(img,1), size(img,2), 3);
        idx = 1;
        idx2 = 0;
        gain = 1;
        cm = lines(size(tr,2));
        cm = brighten(cm,0);
        % cm = flipud(cm);
        value = mean(img,3);
        value = value - min(value(:));
        value = value./(max(value(:))*0.7);
        normalization = zeros(size(value));
        for n = 1:length(space)
            tmp = space{n};
            for p = 1:size(tmp,2)
                if ~ismember(idx,del)
                    idx2 = idx2 + 1;
                    ymin = round(min(a{1,n}(:,1)));
                    ymax = round(max(a{1,n}(:,1)));
                    xmin = round(min(a{1,n}(:,2)));
                    xmax = round(max(a{1,n}(:,2)));
                    ic = tmp(:,p);
                    if sum(ic(:))==0
                        ic(:) = 0.1;
                    else
                        ic = ic-min(ic);
                        ic = ic./max(ic);
                    end
%                     icm = flipud(toimg(ic, (xmax-xmin)+1, (ymax-ymin)+1));
                    icm = toimg(ic, (xmax-xmin)+1, (ymax-ymin)+1);
                    saturation = zeros(size(img,1), size(img,2));
                    saturation(xmin:xmax,ymin:ymax) = icm;
                    valuetmp = value;
                    valuetmp(saturation(:)==0) = 0;
                    tim = zeros(size(value,1), size(value,2), 3);
                    hh = rgb2hsv(cm(idx2,:));
                    tim(:,:,1) = hh(1);
                    tim(:,:,2) = saturation;
                    tim(:,:,3) = valuetmp.*gain;
                    normalization = normalization + logical(valuetmp);
                    tim = hsv2rgb(tim);
                    timf = timf + tim;
                end
                idx = idx + 1;
            end
        end
        valuetmp = value;
        valuetmp(timf(:,:,1)~=0) = 0;
        normalization = normalization + logical(valuetmp);
        
        timf = timf + repmat(valuetmp,1, 1, 3);
        timf = timf./repmat(normalization,1, 1, 3);
        figure
        subplot(4,1,1)
        imshow2(timf)
        f1 = subplot(4,1,2);
        hold on
        for n = 1:size(tr,2)
%             if mean(tr(:,n))>0
%                 tr(:,n) = tr(:,n).*-1;
%             end
            plot(xAxImg(500:size(tr,1)+499),tr(:,n)+n*max(tr(:)),'Color', cm(n,:))
        end
        f2 = subplot(4,1,3);
        plot(xAxVNR,stdV)
        hold on
        plot(xAxVNR(bursts),stdV(bursts),'.');
        f3 = subplot(4,1,4);
        plot(xAxVNR(bursts(2:end)),instFreq,'o');
        ylabel('inst. frequency / Hz')
        xlabel('time/s')
        axis tight
        linkaxes([f1 f2 f3], 'x')
        saveas(gcf,strrep(files(k).name,'meanFrame.tif','pca_ica.fig'))
    end
end