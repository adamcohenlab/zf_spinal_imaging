function [av, pk2, pkXc, f, indivInstFreq, phb, burstIdx] = calcVNRphase(burstImg, stdV, tr, instFreq)
ba = burstImg-499;
bi = 1:length(ba);
noBurst = ba<1 | ba>size(tr,1);
ba(noBurst) = [];
bi(noBurst) = [];
instFreq(noBurst) = [];
phb = [];
idx = 1;
indPeak = [];
indivInstFreq = [];
burstIdx = [];
% smTr = smooth(tr-medfilt1(tr,201),10,'moving');
smTr = tr-medfilt1(tr,201);
for n = 1:length(ba)-1
    if instFreq(n+1)>15 && instFreq(n+1)<40
        interval = ba(n):ba(n+1);
        indivInstFreq(idx) = instFreq(n+1);
        newInterval = linspace(interval(1),interval(end),100);
        tmpTr = interp1(interval,smTr(interval,:),newInterval,'pchip');
        [~,indPeak(idx)] = max(tmpTr);
        phb(:,idx)  = tmpTr;
        burstIdx(idx) = bi(n);
        idx = idx + 1;
    end
end


av = squeeze(mean(phb,2));


% pk2 = angle(mean(exp(1i*pi*(indPeak./50)),1))./pi;
% 
% pk2 = mod((pk2+2),2)*50;
pk2 = indPeak;

%hack to get sampling rate
srate = size(stdV,1)./25;

interpStdV = interp1(1/srate:1/srate:25,stdV,1/996.4:1/996.4:25);

smTr = smooth(tr,10,'moving');

xc = xcorr(interpStdV(500:end),smTr,50);

[~,pkXc] = max(xc);
pkXc = pkXc-50;


[pw,px] = pspectrum(interpStdV,996.4);
[~,fIdx] = max(pw(px>15));
f  = px(px>15);
f = f(fIdx);
end