function spikes = findSpikes(tr)
spikes = cell(size(tr,2),1);
for n = 1:size(tr,2)
    cellN = n;
    %     tt = tr(:,cellN) - medfilt1(tr(:,cellN),21);
    if sum(isnan(tr(:,cellN)))==0
        tt = highpass(tr(:,cellN), 0.04);
    else
        tt = tr(:,cellN);
    end
    % figure
    % plot(tt)
    [N,edges] = histcounts(tt,200);   % histogram of std trace
    %     [~, maxIdx] = max(N);                       % histogram peak
    %     noise = [N(1:maxIdx) fliplr(N(1:maxIdx-1))];% mirror left side of histogram
    %     fobj = fit((1:length(noise))', noise', 'gauss1');   % fit gausian
    fobj = fit((1:length(N))', N', 'gauss1');   % fit gausian
    EdgeIdx = round(4*fobj.c1/sqrt(2) +fobj.b1);
    if EdgeIdx<=size(edges,2)
        thr = edges(round(4*fobj.c1/sqrt(2) +fobj.b1));
    else
        thr = edges(end);
    end% peak threshold is 5*std of gaussian fit
    [peaks,peakInd] = findpeaks(tt,'MinPeakHeight',thr, 'MinPeakDist', 5);
    snr = mean(peaks)./std(tt(~ismember(1:length(tt),peakInd)));
%     figure(444); clf
%     plot(tr(:,cellN))
%     legend(['snr = ' num2str(snr)])
%     hold on
%     plot(peakInd,tr(peakInd,cellN),'.r')
%     inp = input('ok? ','s');
%     if strcmp(inp,'y')
        spikes{cellN} = peakInd;
%     else
%         spikes{cellN} = [];
%     end
end
end