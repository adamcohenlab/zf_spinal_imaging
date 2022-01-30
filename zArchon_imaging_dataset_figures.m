% load data
load('..\Data\zArchonData20210305.mat')
%% 
nFish = length(unique([data([data.brainRegion]=='sc').fishId]))
nCells = length([data([data.brainRegion]=='sc' & [data.trial]==1).rois])
nCpT = [];
id = 1;
for n = unique([data([data.brainRegion]=='sc' & [data.trial]==1).stdV])
    nCpT(id) = length([data([data.brainRegion]=='sc' & [data.trial]==1 & [data.stdV] == n).rois]);
    id = id + 1;
end
meanCells = mean(nCpT)
stdCells = std(nCpT)

srImg = 996.3;
xAxImg = 1/srImg:1/srImg:25000/srImg;
srate = 50000;
xAxVNR = 1/srate:1/srate:length(stdVAll{1})/srate;

selector = false(size(data));
for n = 1:length(data)
    selector(n) = data(n).trial==1 &&...
        data(n).brainRegion=='sc' &&...
        ~isempty(data(n).av) &&...
        ~isempty(data(n).xyzpos) &&...
        ~isnan(data(n).xyzpos(:,3));
end
tr = [data(selector).trace];
masterRef = cat(1,data(selector).stdV);
ROIs = cat(1,data(selector).rois);
bsi = {data(selector).boutStartIdx};
bei = {data(selector).boutEndIdx};
bImg = {data(selector).burstImg};
av = cat(2,data(selector).av);
amp = max(av)-min(av);
ids = find(selector);



% Fig 1C example recording
f = figure;
f.Renderer = 'painters';
tiledlayout(15,1, 'Tilespacing', 'none', 'Padding', 'none')
cm  = cmocean('-matter', 100);
idx = 0;
cells = [1 2 3 5 6 7 8 9 10 11 12 13 14];
[~, si] = sort(amp(cells));
t1 = nexttile(1, [13 1]);
cls = find(selector);
cls(cells)
speedCorrect = data(cls(1)).gainAdjust*data(cls(1)).omrGain;
hold on
for n = 1:length(si)
    if amp(cells(si(n)))>=.5
        plot(xAxImg(500:size(tr,1)+499), tr(:,cells(si(n)))+idx, 'Color', cm(70,:))
    else
        plot(xAxImg(500:size(tr,1)+499), tr(:,cells(si(n)))+idx, 'Color', cm(10,:))
    end
    idx = idx + 6.2;
end
t2 = nexttile(14);
plot(xAxVNR, stdVAll{1}, 'k')
t3 = nexttile(15);
[dat,~] = import2pdaq('helper_data\', 'l1-tr1.A2','a');
plot(xAxVNR, medfilt1(dat, 401).*speedCorrect.*1.3.*20, 'k')
linkaxes([t1 t2 t3], 'x')
xlim([4 6.5])

% Fig 1F phase vs position
imgMat = loadtiff('helper_data\rec00005_meanFrame.tif');
V2as = ids(amp>=0.5);
stackId = cat(1,data(V2as).fishId).*10  + cat(1,data(V2as).xypos);
recId = cat(1,data(V2as).stdV);
uniqueRecId = unique(recId);
uniqueId = unique(stackId);
xPos = cat(1, data(V2as).xyzpos);
xPos = xPos(:,1);
xScale = (100/(180/25)/6.5); % f_tubelens / (f_tubOlympus/M_objective)/pixel_size)

%use individual bouts
tmp = [];
tmp2 = [];
phaseInd = [];
f = [];
frm_std = [];
for n = 1:length(data)
    
    %correct for shift due to camera row clock
    pk2Corr = data(n).pk2;
    indPeriod = data(n).selInstFreq;
    indPeriod = 1./indPeriod;
    pkTime = pk2Corr.*indPeriod./100;
    % distance from center * row delay in seconds
    if ~isempty(data(n).xyzpos)
        pkTime = pkTime + abs(data(n).xyzpos(2)-(size(imgMat,1)/2 - 1))*9.74436e-6;
        pk2Corr = pkTime.*100./indPeriod;
        phaseMean = circ_mean(2.*pi.*pk2Corr'./100);
        phaseMean = 100*phaseMean/(2*pi);
    else
        phaseMean = [];
    end
    
        [~, frm] = max(data(n).av);
        if ~isempty(exp(1i.*pi.*(frm./50))) && ~isempty(data(n).xyzpos)
            tmp(n) = exp(1i.*pi.*(phaseMean./50));
%             tmp(n) = exp(1i.*pi.*(frm./50));
            pos = data(n).pk2;
            frm_std(n) = (circ_std(2.*pi.*pos'./100))/sqrt(numel(pos));
        else
            tmp(n) = NaN;
        end


    
    
    tmp2{n} = exp(1i.*pi.*(pk2Corr./50));
    phaseInd{n} = exp(1i.*(angle(tmp2{n}).*2));
    instf = data(n).instFreq; 
    ba = data(n).burstImg-499;
    noBurst = ba<1 | ba>size(data(n).trace,1);
    instf(noBurst) = [];
    f{n} = instf(instf>15 & instf<40);
end
phaseInd = phaseInd(V2as);
f = f(V2as);
frmStd = frm_std(V2as);
phaseFrames = cat(1,tmp(V2as));
phase = exp(1i.*(angle(phaseFrames).*2))';
phase2 = phaseFrames';

avSel = cat(2, data(V2as).av);

% figure
% hold on
plts = ceil(sqrt(length(uniqueId)));
allPhase = [];
allPhase2 = [];
avShift = [];
xP = [];
xpi = [];
recId2 = [];
stdAll = [];
for n = 1:length(uniqueId)
    if sum(stackId==uniqueId(n))>1
        tmpX = xPos(stackId==uniqueId(n));
        tmpStd = frmStd(stackId==uniqueId(n));
        tmpPhase = phase(stackId==uniqueId(n));
        tmpPhase2 = phase2(stackId==uniqueId(n));
        avPhase = mean(tmpPhase);
        
        avShift = cat(2, avShift, circshift(avSel(:,stackId==uniqueId(n)),...
            round((angle(avPhase)/pi)*50)));
        phaseRot = exp(1i.*(angle(tmpPhase)-angle(avPhase)));
        phase2Rot = exp(1i.*(angle(tmpPhase2)-angle(avPhase)/2-.5*pi));
        
        allPhase = cat(1,allPhase, phaseRot);
        allPhase2 = cat(1, allPhase2, phase2Rot);
        recId2 = cat(1, recId2, recId(stackId==uniqueId(n)));
%         subplot(plts,plts,n)
%         scatter(tmpX, angle(tmpPhase)./pi,...
%             10, 'filled')
%         hold on
%         line([tmpX tmpX]', [(angle(tmpPhase)'-tmpStd./2)./pi; (angle(tmpPhase)'+tmpStd./2)./pi], 'Color', 'k')
%         title(num2str(angle(avPhase)/pi))
%         xlim([0 1000])
%         ylim([-1 1])
        xP = cat(1,xP, tmpX);
        stdAll = cat(2,stdAll, tmpStd);
        xpi = cat(1,xpi, true(size(tmpX)));
    else
        xpi = cat(1,xpi, false);
    end
end




%Fig 1H
figure
hold on
xVal = xP./xScale;
yVal = angle(allPhase2)./pi;
xVal(yVal>0) = xVal(yVal>0)+500;

ft = fittype(@(a,b1,b2,x) (a*x + b1).*(x<500) + (a*x + b2).*(x>=500));
[fobj, gof, out] = fit(xVal, yVal, ft);

alpha = 0.95
ci = confint(fobj, alpha) 
t = tinv((1+alpha)/2, gof.dfe); 
se = (ci(2,:)-ci(1,:)) ./ (2*t) % Standard Error




xVal = xP./xScale;
scatter(xVal, yVal, 10, 'b', 'filled')
plot(0:450, fobj(0:450), 'r')
plot(0:450, fobj(500:500+450), 'r')
legend('neurons', ['slope = ' num2str(fobj.a)])
xlabel('distance (\mum)')
ylabel('phase (\pi)')


wavelength_in_mm = 2/(fobj.a*1000)

% F1 phase map

%Fig 1G
figure
yzPos = cat(1, data(ids).yzproj);
xPos = cat(1, data(ids).xyzpos);
xPos = xPos(:,1);
bgColor = [1 1 1];
set(gcf, 'Color', bgColor)
subplot(1,2,1)
scatter(yzPos(:,1), yzPos(:,2).*-1 +1, 30, [0.8 0.8 0.8],...
    'filled', 'MarkerFaceAlpha', 1)
hold on
yzPos = cat(1, data(V2as).yzproj);
yzPos = yzPos(xpi==1,:);
xPos = xP;
phaseIdx = angle(allPhase2)./pi;
phaseIdx = (phaseIdx.*49.5 + 0.5) +50;
phaseIdx = round(phaseIdx);
bgColor = [1 1 1];
plotAlpha = 1;
C = cmocean('phase', 100);
subplot(1,2,1)
scatter(yzPos(:,1), yzPos(:,2).*-1 +1, 30, C(phaseIdx,:),...
    'filled', 'MarkerFaceAlpha', plotAlpha)
hold on
plot(0.5.*sin(linspace(0,2*pi,100)) + 0.5,...
    0.5.*cos(linspace(0,2*pi,100)) + 0.5,'k')
set(gca, 'DataAspectRatio', [1,1,1])

set(gca, 'Color', bgColor)
set(gca, 'XColor', 'k')
set(gca, 'YColor', 'k')
ylim([0 1])
zlim([0 1])
grid('off')
subplot(1,2,2)
scatter3(xPos, yzPos(:,1), yzPos(:,2).*-1 +1, 30, C(phaseIdx,:),...
    'filled', 'MarkerFaceAlpha', plotAlpha)
set(gca, 'DataAspectRatio', [1024,1,1])
set(gca, 'Color', bgColor)
set(gca, 'XColor', 'k')
set(gca, 'YColor', 'k')
ylim([0 1])
view(0,90)
grid('off')
axis off

clear cbar;
cbar(1,:,1) = ones(100,1).*C(:,1);
cbar(1,:,2) = ones(100,1).*C(:,2);
cbar(1,:,3) = ones(100,1).*C(:,3);
cbar = repmat(cbar,5, 1, 1);

axes('Position', [0.1 0.1 .35 .1])
imshow2(cbar)
text(0,-5,num2str(0,2), 'HorizontalAlignment', 'center', 'Color', 'k')
text(100,-5,'\pi', 'HorizontalAlignment', 'center', 'Color', 'k')
set(gcf, 'Renderer', 'painters')
% phase example
exampleSelector = [data.fishId]==4 & [data.zPos]==2 & [data.trial]==1;
exampleIdx = find(exampleSelector);
exampleTraces = cat(2, data(exampleSelector).trace);
examplePhases = angle(tmp(exampleSelector));
examplePhases = mod(examplePhases-nanmean(examplePhases), 2*pi);
examplePhases = ceil(99.*(examplePhases)./(2.*pi));
% examplePhases = ceil(99.*(examplePhases + pi)./(2.*pi));

%Fig 1F
figure
tiledlayout(2,1, 'Padding', 'none')
nexttile(2)
hold on
C = cmocean('phase', 100);
idx = 0;
[~, si] = sort(examplePhases);
% si = 1:31;
for n = 1:size(exampleTraces, 2)
    if ~isnan(exampleTraces(1,si(n)))
        tmpAmp = max(data(exampleIdx(si(n))).av) - min(data(exampleIdx(si(n))).av);
        if tmpAmp>=.5
            tmpColor = C(examplePhases(si(n)),:);
        else
            tmpColor = [0 0 0];
        end
        plot(xAxImg(1:length(exampleTraces(:,1))), exampleTraces(:,si(n))+idx.*4, 'Color', tmpColor)
        plot(xAxImg(data(exampleIdx(si(n))).spikes), ones(size(data(exampleIdx(si(n))).spikes))+idx.*4, '.k')
        idx = idx + 1;
    end
end
xlabel('time (s)')

load('helper_data\rec00005_ica_pca_selected.mat', 'space');
% imgMat = loadtiff('X:\Lab\Labmembers\Urs Boehm\Data\7-zArchon\20190228\l2\rec00005_meanFrame.tif');
timf = zeros(size(imgMat,1), size(imgMat,2), 3);
for n = 1:size(exampleTraces, 2)
    if ~isnan(exampleTraces(1,n))
        
        roi = data(exampleIdx(n)).rois{1};
        ymin = round(min(roi(:,1)));
        ymax = min([round(max(roi(:,1))) size(imgMat,2)]);
        xmin = round(min(roi(:,2)));
        xmax = min([round(max(roi(:,2))) size(imgMat,1)]);
        ic = space{n};
        ic = ic-min(ic);
        ic = ic./max(ic);
        icm = flipud(toimg(ic, (xmax-xmin)+1, (ymax-ymin)+1));
        saturation = zeros(size(imgMat,1), size(imgMat,2));
        saturation(xmin:xmax,ymin:ymax) = icm;
        tmpAmp = max(data(exampleIdx(n)).av) - min(data(exampleIdx(n)).av);
         if tmpAmp>=.5
            tmpColor = C(examplePhases(n),:);
        else
            tmpColor = [1 1 1];
         end
        tim = zeros(size(imgMat,1), size(imgMat,2), 3);
        tim(:,:,1) = saturation.*tmpColor(1);
        tim(:,:,2) = saturation.*tmpColor(2);
        tim(:,:,3) = saturation.*tmpColor(3);
        timf = timf + tim;
        
    end
end

nexttile(1)
imagesc([0 size(timf,2)/xScale], [1 size(timf,1)/xScale], timf)
daspect([1 1 1])
hold on
xlabel('space (\mum)')
for n = 1:size(exampleTraces, 2)
    if ~isnan(exampleTraces(1,si(n)))
        roi = data(exampleIdx(si(n))).rois{1};
        [x, y] = centroid(polyshape(roi));
        text(x/xScale, y/xScale, num2str(n), 'Color', 'w')
    end
end

% Fig 1D
sID = [data.fishId]==4 & [data.zPos]==2 & [data.trial]==1;
exampleId = find(sID(selector));

C = cmocean('-matter', 100);
srImg = 996.3;
xAxImg = 1/srImg:1/srImg:25000/srImg;
srate = 50000;
xAxVNR = 1/srate:1/srate:length(stdVAll{1})/srate;
exampleTrace = find(selector);
f1 = figure;
f1.Renderer = 'painters';
s1 = axes('OuterPosition', [0, 0.6, 0.8, 0.3]);
hold on
plot(xAxImg(500:size(tr,1)+499), tr(:,exampleId(21)), 'k')

ylim([-1.5 4])
s2 = axes('OuterPosition', [0, 0.3, 0.8, 0.3]);
hold on
plot(xAxImg(500:size(tr,1)+499), tr(:,exampleId(23)), 'k')
ylim([-1.5 4])
s3 = axes('OuterPosition', [0, 0, 0.8, 0.3]);
plot(xAxVNR, stdVAll{data(exampleTrace(exampleId(21))).stdV}, 'k')
hold on
plot(xAxVNR, stdVAll{data(exampleTrace(exampleId(21))).stdV}, 'k')
plot(xAxVNR(data(exampleTrace(exampleId(21))).burstVNR), stdVAll{data(exampleTrace(exampleId(21))).stdV}(data(exampleTrace(exampleId(21))).burstVNR), 'or')
linkaxes([s1 s2 s3], 'x')
xlim([6.55 6.95])
ylim([0.00 0.02])

ampIdx = amp;
ampIdx = (ampIdx - min(ampIdx))/max(ampIdx - min(ampIdx));
ampIdx = ampIdx.*99 +1;
ampIdx = round(ampIdx);

axes('OuterPosition', [0.8, 0.6, 0.2 0.3])
hold on
p1 = plot(xAxImg(1:100), data(exampleTrace(exampleId(21))).indCycle, 'Color', [0 0 0 0.1]);
plot(xAxImg(1:100), av(:,exampleId(21)), 'Color', C(ampIdx(exampleId(21)),:), 'LineWidth', 2)
ylim([-1.5 4])
axes('OuterPosition', [0.8, 0.3, 0.2 0.3])
hold on
plot(xAxImg(1:100), data(exampleTrace(exampleId(23))).indCycle, 'Color', [0 0 0 0.1])
plot(xAxImg(1:100), av(:,exampleId(23)), 'Color', C(ampIdx(exampleId(23)),:), 'LineWidth', 2)
ylim([-1.5 4])
axes('OuterPosition', [0.8, 0, 0.2 0.3])
hold on
bursts = data(exampleTrace(exampleId(21))).burstVNR;
avVNRburst = zeros(100, length(bursts)-1);
tmpVNR = stdVAll{data(exampleTrace(exampleId(21))).stdV};
for n = 1:length(bursts)-1
    tmpBurst = tmpVNR(bursts(n):bursts(n+1));
    avVNRburst(:,n) = interp1(1/length(tmpBurst):1/length(tmpBurst):1, tmpBurst, 1/100:1/100:1);
end
plot(avVNRburst, 'Color', [0 0 0 0.1])
plot(mean(avVNRburst,2), 'Color', [0 0.4 1], 'LineWidth', 2)
ylim([0.00 0.02])

%Fig 2B
% spike probability over burst cycle in v3s
freqCutoff = [10 40];
selector = false(size(data));
for n = 1:length(data)
    selector(n) = data(n).trial==1 &&...
        data(n).brainRegion=='sc' &&...
        ~isempty(data(n).av) &&...
        ~isempty(data(n).xyzpos) &&...
        ~isnan(data(n).xyzpos(:,3));
end
tr = [data(selector).trace];
av = cat(2,data(selector).av);
amp = max(av)-min(av);
ids = find(selector);
yzPos = cat(1, data(ids).yzproj);
v3 = ids(yzPos(:,2)>0.7);
stackId = cat(1,data(v3).fishId).*10  + cat(1,data(v3).xypos);
uniqueId = unique(stackId);
xPos = cat(1, data(v3).xyzpos);
xPos = xPos(:,1);

boutSpikeRate = [];
interboutSpikeRate = [];
boutVm = [];
interboutVm = [];
boutMod = [];
swimmingAll = [];
for n = 1:length(data)
    boutSpikes = 0;
    totBouts = 0;
    interboutSpikes = 0;
    totInterbouts = 0;
    for trls = n%find([data.cellId]==data(n).cellId)
        ba = data(trls).burstImg-499;
        bs = ba(data(trls).boutStartIdx);
        be = ba(data(trls).boutEndIdx);
        %get rid of bouts that are outside the imaging time
        noBout = bs<1 | be<1 | bs>size(tr,1) | be>size(tr,1);
        bs(noBout) = [];
        be(noBout) = [];
        noBurst = ba<bs(1) | ba>be(end);
        ba(noBurst) = [];
        
        [~,bsi] = ismember(bs,ba);
        [~,bei] = ismember(be,ba);
        swimming = zeros(size(tr(:,1)));
        for k = 1:length(bs)
            swimming(bs(k):be(k)) = 1;
        end
        boutSpikes = boutSpikes + sum(swimming(data(n).spikes));
        totBouts = totBouts + sum(swimming);
        interboutSpikes = interboutSpikes + sum(~swimming(data(n).spikes));
        totInterbouts = totInterbouts + sum(~swimming);
    end
    swimmingAll(:,n) = swimming;
    boutSpikeRate(n) = boutSpikes/totBouts;
    interboutSpikeRate(n) = interboutSpikes/totInterbouts;
end

hcSpks = [];
spksAll = {};
vnrBurstMean = [];
for k = 1:length(data)
    spksSc = [];
    vnrTmp = [];
    idx = 1;
    %trls = k, each trial separate, trls = find([data..., all trials of
    %cell combined
    for trls = k%find([data.cellId]==data(k).cellId)
        ba = data(trls).burstImg-499;
        bv = data(trls).burstVNR;
        noBurst = ba<1 | ba>size(tr,1);
        ba(noBurst) = [];
        bv(noBurst) = [];
        instFreq = data(trls).instFreq;
        instFreq(noBurst) = [];
        tmpSpk = {};
        for n = 1:length(ba)-1
            if instFreq(n+1)>freqCutoff(1) && instFreq(n+1)<freqCutoff(2)
                interval = ba(n):ba(n+1);
                spks = data(trls).spikes;
                spks = (spks(spks>=ba(n) & spks<ba(n+1)));
                tmpSpk{n} = (spks-ba(n))/(ba(n+1)-ba(n));
                spksSc = cat(1, spksSc, (spks-ba(n))/(ba(n+1)-ba(n)));
                idx = idx + 1;
                tmpBurst = stdVAll{data(trls).stdV}(bv(n):bv(n+1));
                vnrTmp = cat(1, vnrTmp, interp1(1/length(tmpBurst):1/length(tmpBurst):1, tmpBurst, 1/100:1/100:1));
            end
        end
    end
    spksAll{k} = tmpSpk;
    hcSpks(:,k) = histcounts(spksSc, linspace(0,1,21))';
    hcSpks(:,k) = hcSpks(:,k)/idx;
    hcSpks(:,k) = hcSpks(:,k)./(1/20);
    vnrBurstMean(:,k) = mean(vnrTmp,1);
end

logSR = (boutSpikeRate-interboutSpikeRate)./(boutSpikeRate+interboutSpikeRate);
logSR(isnan(logSR)) = 0;


sel = v3(ismember(v3, find(logSR>0)));
figure
tiledlayout(2,1)
nexttile
hold on
xa = linspace(0,1,20);
av = mean(hcSpks(:,sel),2);
sem = std(hcSpks(:,sel),[],2)./sqrt(size(hcSpks(:,sel),2));
patch([xa fliplr(xa)], [av+sem; flipud(av-sem)]', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5)
plot(xa,av, 'k', 'Linewidth', 2)
ylim([0 0.6])
ylabel('mean spike probability')
xlabel('burst phase')
nexttile
plot(mean(vnrBurstMean(:,sel),2), 'Color', [0 0.4 1], 'LineWidth', 2)


% Fig2B
figure
hold on
id = 0;
for n = 1:size(spksAll{v3(1)},2)
    x = spksAll{v3(1)}{n};
    if ~isempty(x)
    for k = 1:length(x)
%         line([x(k) x(k)], [id id+1], 'Color', 'k')
        plot([x(k)], [id], '.k')
    end
    id = id + 1;
    end
end

%Fig S3 A
% logSR map
yzPos = cat(1, data(ids).yzproj);
figure
set(gcf, 'Color', bgColor)
hold on
subplot(1,2,1)
scatter(yzPos(:,1), yzPos(:,2).*-1 +1, 30, [0.8 0.8 0.8],...
    'filled', 'MarkerFaceAlpha', 1)
hold on


nonOsc = ids(amp<.5);

yzPos = cat(1, data(nonOsc).yzproj);
xPos = cat(1, data(nonOsc).xyzpos);

nns = isnan(xPos(:,3));
xPos = xPos(:,1);
xPos(nns) = [];
ampIdx = logSR(nonOsc);
% ampIdx(ampIdx>.6) = .6;
% ampIdx(ampIdx<-.6) = -.6;
ampIdx(nns) = [];
maxAbs = max(sqrt(ampIdx.^2));
ampIdx(end+1) = maxAbs;
ampIdx(end+1) = -1*maxAbs;
% ampIdx(isinf(ampIdx)) = max(ampIdx(~isinf(ampIdx)));
ampIdx = (ampIdx - min(ampIdx))/max(ampIdx - min(ampIdx));
ampIdx = ampIdx.*99 +1;
ampIdx = round(ampIdx);
ampIdx(end-1:end) = [];
bgColor = [1 1 1];
C = cmocean('balance', 100);
% C = C.*0.9;

% axes('Position',[0.05 0.05 .25 1])
scatter(yzPos(:,1), yzPos(:,2).*-1 +1, 30, C(ampIdx,:),...
    'filled', 'MarkerFaceAlpha', 1)
hold on
plot(0.5.*sin(linspace(0,2*pi,100)) + 0.5,...
    0.5.*cos(linspace(0,2*pi,100)) + 0.5,'k')
% view(90,0)
set(gca, 'DataAspectRatio', [1,1,1])
set(gca, 'Color', bgColor)
% xlim([0 1])
ylim([0 1])
zlim([0 1])
grid('off')
% axes('Position',[.4 0.05 .65 1])
subplot(1,2,2)
scatter3(xPos, yzPos(:,1), yzPos(:,2).*-1 +1, 30, C(ampIdx,:),...
    'filled', 'MarkerFaceAlpha', 0.5)


set(gca, 'DataAspectRatio', [1024,1,1])
set(gca, 'Color', bgColor)
view(0,0)
grid('off')

clear cbar;
cbar(1,:,1) = ones(100,1).*C(:,1);
cbar(1,:,2) = ones(100,1).*C(:,2);
cbar(1,:,3) = ones(100,1).*C(:,3);
cbar = repmat(cbar,5, 1, 1);

axes('Position', [0.1 0.1 .35 .1])
imshow2(cbar)
text(0,-5,['<' num2str(-1*maxAbs,2)], 'HorizontalAlignment', 'center')
text(100,-5,['>' num2str(maxAbs,2)], 'HorizontalAlignment', 'center')


%Fig. S3B
nonOscTraces = cat(2,data(nonOsc).trace);
[~, sri] = sort(logSR(nonOsc));
t1 = find((-1.*yzPos(:,2)+1)>=0.887205 & (-1.*yzPos(:,2)+1)<0.887206);
t2 = find((-1.*yzPos(:,2)+1)>=0.824 & (-1.*yzPos(:,2)+1)<0.824119);
t3 = find((-1.*yzPos(:,2)+1)>=0.748682 & (-1.*yzPos(:,2)+1)<0.748683);
t4 = find((-1.*yzPos(:,2)+1)>=0.635136 & (-1.*yzPos(:,2)+1)<0.635137);
% t3 = find((-1.*yzPos(:,2)+1)>=0.7302 & (-1.*yzPos(:,2)+1)<0.730222);
figure
tiledlayout(2,1)
nexttile
hold on
plot(xAxImg(501:500+size(nonOscTraces,1)),nonOscTraces(:,[t1]), 'Color', C(ampIdx(t1),:))
plot(xAxImg(501:500+size(nonOscTraces,1)),nonOscTraces(:,[t2])+6, 'Color', C(ampIdx(t2),:))
a1 = gca;
nexttile
plot(xAxVNR, stdVAll{data(nonOsc(t1)).stdV}, 'k')
a2 = gca;
linkaxes([a1 a2], 'x')

%Fig 1E
% pie charts dorsal vs ventral
dorsalNeurons = ~ismember(ids, v3);
ventralNeurons = ismember(ids, v3);
nDorsal = sum(dorsalNeurons);
nVentral = sum(ventralNeurons);
% oscilating swim-active nonswim-active
dorsal = [sum(amp(dorsalNeurons)>.5),...
    sum(amp(dorsalNeurons)<=.5 & logSR(ids(dorsalNeurons))>0),...
    sum(amp(dorsalNeurons)<=.5 & logSR(ids(dorsalNeurons))<=0)];

ventral = [sum(amp(ventralNeurons)>.5),...
    sum(amp(ventralNeurons)<=.5 & logSR(ids(ventralNeurons))>0),...
    sum(amp(ventralNeurons)<=.5 & logSR(ids(ventralNeurons))<=0)];

figure
tiledlayout(2,1, 'TileSpacing','compact')
nexttile(1)
pd = pie(dorsal, {num2str(dorsal(1)), num2str(dorsal(2)), num2str(dorsal(3))});
nexttile(2)
pv = pie(ventral, {num2str(ventral(1)), num2str(ventral(2)), num2str(ventral(3))});
lgnd = {'oscilating', 'non-oscilating swim-active', 'non-oscilating non-swimactive'};
lgd = legend(lgnd);
lgd.Layout.Tile = 'east';

cm1 = cmocean('-matter', 100);
cm2 = cmocean('balance', 100);
pd(1).FaceColor = cm1(98,:);
pd(1).EdgeColor = 'none';
pd(3).FaceColor = cm2(85,:);
pd(3).EdgeColor = 'none';
pd(5).FaceColor = cm2(15,:);
pd(5).EdgeColor = 'none';
pv(1).FaceColor = cm1(98,:);
pv(1).EdgeColor = 'none';
pv(3).FaceColor = cm2(85,:);
pv(3).EdgeColor = 'none';
pv(5).FaceColor = cm2(15,:);
pv(5).EdgeColor = 'none';

% ypos vs bout amp scatter
%Fig 1E
colors = zeros(size(amp));
colors(amp>.5) = 1;
colors(amp<=.5 & logSR(ids)>0) = 2;
colors(amp<=.5 & logSR(ids)<=0) = 3;
cm = [cm1(98,:); cm2(85,:); cm2(15,:)];
figure
yzPos = cat(1, data(ids).yzproj);
scatter(amp, 1-yzPos(:,2), 'filled', 'SizeData', 5, 'CData', cm(colors,:))
ylabel('d-v position')
xlabel('VNR triggered average amplitude')
%% setup longTable for amp vs SR scatter
[cellIds, ia, ic] = unique([data.cellId]);
zPos = nan(size(cellIds));
avAmp = nan(size(cellIds));
for n = 1:length(cellIds)
    if ~isempty(data(ia(n)).yzproj)
        zPos(n) = data(ia(n)).yzproj(2);
        avAmp(n) = max(data(ia(n)).av)-min(data(ia(n)).av);
    end
end
sw = 0; 
for n = 1:length(data)
    sw = sw+length(data(n).boutStartIdx);
end
longTable = table('Size', [sw, 26], 'variableTypes', {'double', 'double', 'double',...
    'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double',...
    'double', 'double', 'double', 'double', 'double', 'double', 'cell', 'cell', 'cell', 'cell', 'logical'},...
    'variableNames', {'NSpikes', 'dataN', 'BoutDuration', 'InterboutDuration', 'BoutAmp', 'Speed', 'SpeedRaw', 'DistanceRaw', 'DistanceReal', 'SpeedCorrected', 'SpeedRawTot', 'Ypos',...
    'CycleAmp', 'CellId', 'FishId', 'FOV', 'BoutId', 'Trial', 'Freq', 'recId', 'BoutAv',...
    'SpikesEnd', 'Bursts', 'Spikes', 'BurstsEnd', 'UseFlag'});
id = 1;
currVnr = 0;
currFov = 1;
fovK = 0;
for n = 1:length(data)
    if currVnr ~= data(n).stdV
        tmpVnr = stdVAll{data(n).stdV};
        currVnr = data(n).stdV;
        sr = length(tmpVnr)/25;
        xax = 1/sr:1/sr:25;
        tmpVnr = timeseries(tmpVnr, xax);
        tmpVnr = idealfilter(tmpVnr, [0 5], 'pass');
        tmpVnr = tmpVnr.Data;
        tmpVnr = tmpVnr-min(tmpVnr);
        tmpFlow = gratingFlow{data(n).stdV};
    end
    if n>1 && (data(n).xypos~=data(n-1).xypos || data(n).fishId~=data(n-1).fishId || data(n).zPos~=data(n-1).zPos)
        currFov = currFov + 1;
    end
    tmpTr = data(n).trace;
    selBursts = find(data(n).instFreq>15 & data(n).instFreq<40 & data(n).burstImg'>499 & data(n).burstImg'<=length(tmpTr));
    selBursts = selBursts(2:end);
    
    if n>1 && data(n).stdV~=data(n-1).stdV
        fovK = longTable.BoutId(id-1);
    end
    for k = 1:length(data(n).boutStartIdx)
      
        longTable.BoutId(id) = k+fovK;
        %data(n).burstImg is untrimmed, data(n).spikes is trimmed
        longTable.NSpikes(id) = sum((data(n).spikes+500)>= data(n).burstImg(data(n).boutStartIdx(k)) & ...
            (data(n).spikes+500)<= data(n).burstImg(data(n).boutEndIdx(k)));
        %don't use any bouts that are outside of the trimmed region
        if data(n).burstImg(data(n).boutStartIdx(k))>499 && data(n).burstImg(data(n).boutEndIdx(k))<=length(data(n).trace)
            longTable.UseFlag(id) = true;
        else
            longTable.UseFlag(id) = false;
        end
        
        bouts = find(selBursts>=data(n).boutStartIdx(k) & selBursts<=data(n).boutEndIdx(k));
        if ~isempty(data(n).indCycle)
            tmpAv = mean(data(n).indCycle(:,bouts),2);
        else
            tmpAv = NaN;
        end
        longTable.BoutAv(id) = max(tmpAv)-min(tmpAv);
    
        
        spikePos = data(n).spikes((data(n).spikes+500)>= data(n).burstImg(data(n).boutStartIdx(k)) & ...
            (data(n).spikes+500)<= data(n).burstImg(data(n).boutEndIdx(k))) + 500;
        if ~isempty(spikePos)
            spikePos = data(n).spikes + 500 - data(n).burstImg(data(n).boutStartIdx(k));
            spikePosEnd = data(n).spikes + 500 - data(n).burstImg(data(n).boutEndIdx(k));
            longTable.Spikes(id) = {spikePos};
            longTable.SpikesEnd(id) = {spikePosEnd};
        else
            longTable.Spikes(id) = {NaN};
        end
        
        longTable.Bursts(id) = {data(n).burstImg - data(n).burstImg(data(n).boutStartIdx(k))};
        longTable.BurstsEnd(id) = {data(n).burstImg - data(n).burstImg(data(n).boutEndIdx(k))};
        longTable.BoutDuration(id) = (data(n).burstVNR(data(n).boutEndIdx(k)) -...
        data(n).burstVNR(data(n).boutStartIdx(k)))./(length(tmpVnr)./25);
        if k < length(data(n).boutStartIdx)
            longTable.InterboutDuration(id) = (data(n).burstVNR(data(n).boutStartIdx(k+1)) -...
            data(n).burstVNR(data(n).boutEndIdx(k)))./(length(tmpVnr)./25);
        else
            longTable.InterboutDuration(id) = NaN;
        end
        longTable.BoutAmp(id) = mean(tmpVnr(data(n).burstVNR(data(n).boutStartIdx(k)):...
            data(n).burstVNR(data(n).boutEndIdx(k))));
        
%         longTable.BoutAmp(id) = mean(tmpVnr(data(n).burstVNR(data(n).boutStartIdx(k):data(n).boutEndIdx(k))));
        
        longTable.Speed(id) = mean(tmpFlow(data(n).burstVNR(data(n).boutStartIdx(k)):...
            data(n).burstVNR(data(n).boutEndIdx(k))));
        longTable.DistanceRaw(id) = sum(tmpFlow(data(n).burstVNR(data(n).boutStartIdx(k)):...
            data(n).burstVNR(data(n).boutEndIdx(k)))); 
        longTable.DistanceReal(id) = sum(tmpFlow(data(n).burstVNR(data(n).boutStartIdx(k)):...
            data(n).burstVNR(data(n).boutEndIdx(k))))*data(n).gainAdjust*data(n).omrGain.*1.3.*20./...
            (length(tmpVnr)./25); %mm (not raw dist anymore);
        longTable.SpeedRaw(id) = longTable.Speed(id);
        if k + 1 <= length(data(n).boutStartIdx)
            endpoint = data(n).burstVNR(data(n).boutStartIdx(k+1));
        else
            endpoint = length(tmpFlow);
        end
        longTable.SpeedRawTot(id) = mean(tmpFlow(data(n).burstVNR(data(n).boutStartIdx(k)):endpoint));
        longTable.SpeedCorrected(id) = longTable.Speed(id)*data(n).gainAdjust;
        longTable.Speed(id) = ((longTable.Speed(id)*data(n).gainAdjust*data(n).omrGain)-data(n).omrSpeed).*1.3.*20; %mm/s
        
        longTable.dataN(id) = n;
        longTable.FOV(id) = currFov;
        if ~isempty(data(n).yzproj) && ~isempty(data(n).av)
            longTable.Ypos(id) = data(n).yzproj(2);
            longTable.CycleAmp(id) = max(data(n).av)-min(data(n).av);
            longTable.CellId(id) = data(n).cellId;
            longTable.FishId(id) = data(n).fishId;
            longTable.Trial(id) = data(n).trial;
            longTable.Freq(id) = mean(data(n).instFreq(data(n).boutStartIdx(k)+1:...
                data(n).boutEndIdx(k)));
            longTable.recId(id) = data(n).stdV;
        else
            longTable.Ypos(id) = NaN;
            longTable.CycleAmp(id) = NaN;
            longTable.CellId(id) = NaN;
            longTable.FishId(id) = NaN;
            longTable.Trial(id) = NaN;
        end
        id = id + 1;
        if mod(id, 100)==0
            disp(num2str(id))
        end
    end
end
%%
[selectedBouts, ia, ~] = unique(longTable.BoutId);
selector = longTable.Freq>freqCutoff(1) &...
    longTable.Freq<freqCutoff(2) & longTable.UseFlag & longTable.SpeedRaw >0.03 &...
    (ismember(1:height(longTable), ia))';

f_av = mean(longTable.Freq(selector))
f_std = std(longTable.Freq(selector))
dur_av = mean(longTable.BoutDuration(selector))
dur_std = std(longTable.BoutDuration(selector))
intdur_av = nanmean(longTable.InterboutDuration(selector))
intdur_std = nanstd(longTable.InterboutDuration(selector))

%stats about number of cells in the dataset
selector = longTable.Freq>freqCutoff(1) &...
    longTable.Freq<freqCutoff(2) & longTable.UseFlag & longTable.SpeedRaw >0.03;

size(unique(longTable.CellId(selector)))
size(unique(longTable.FishId(selector)))
fovs = (unique(longTable.FOV(selector)));
nC = zeros(length(fovs),1);
for n = 1:length(fovs)
    nC(n) = length(unique(longTable.CellId(selector & longTable.FOV==fovs(n))));
end
mean(nC)
std(nC)
max(nC)

%stats about V2as used in calcualting phase
selector = longTable.Freq>freqCutoff(1) &...
    longTable.Freq<freqCutoff(2) & longTable.UseFlag & longTable.SpeedRaw >0.03 &...
    (ismember(1:height(longTable), ia))' & ismember(longTable.CellId, [data(V2as).cellId]);
f_av = mean(longTable.Freq(selector))


selector = ismember(longTable.CellId, selectedCells) & longTable.Freq>freqCutoff(1) &...
    longTable.Freq<freqCutoff(2) & longTable.UseFlag & longTable.SpeedRaw >0.03;
selectedTrials = longTable.recId(selector);
selectedTrials = unique(selectedTrials);

%% speedRaw vs SR per bout (fig. 4B plot1)
steps = 0:.7:1.5;
allNrmX = [];
allNrmY = [];
finalId = [];
useCell = [];
selectedCells = [data(v3(ismember(v3, find(logSR>0)))).cellId];
selector = ismember(longTable.CellId, selectedCells) & longTable.Freq>freqCutoff(1) &...
    longTable.Freq<freqCutoff(2) & longTable.UseFlag & longTable.SpeedRaw >0.03;
selectedFovs = unique(longTable.FOV(selector));
selectedBouts = unique(longTable.BoutId(selector));

indCellVals = zeros(length(steps)-1,length(selectedCells));
indCellValsSum = zeros(2,length(selectedFovs));

indFovVals = zeros(2, length(selectedFovs));

yVals = longTable.SpeedRaw(selector);

xVals = double((longTable.NSpikes(selector)./longTable.BoutDuration(selector))>0);

allCellId = longTable.CellId(selector);
allFishId = longTable.FishId(selector);

for n = 1:length(selectedFovs)
    boutsInFov = longTable.BoutId(selector);
    boutsInFov = boutsInFov(longTable.FOV(selector)==selectedFovs(n));
    

    tmpX = xVals(longTable.FOV(selector)==selectedFovs(n));
    tmpY = yVals(longTable.FOV(selector)==selectedFovs(n));
    

    tmpbin = [];

    nrmX = tmpX;

    nrmY = tmpY./mean(tmpY);
% nrmY = tmpY;
    meanX = zeros(size(unique(boutsInFov)));
    meanY = zeros(size(unique(boutsInFov)));
    bi = 1;
    for k = unique(boutsInFov)'
        meanX(bi) = mean(nrmX(boutsInFov==k));
        meanY(bi) = mean(nrmY(boutsInFov==k));
        bi = bi + 1;
    end
        allNrmX = [allNrmX; meanX];
        allNrmY = [allNrmY; meanY];
        
        for k = 1:2
            indCellValsSum(k,n) = sum(meanX>=(k-1.5) & meanX<(k-0.5));
            if sum(meanX>=k-1.5 & meanX<k-0.5)>1
                indFovVals(k,n) = mean(meanY(meanX>=k-1.5 & meanX<k-0.5));
            else
                indFovVals(k,n) = NaN;

            end
           
        end
end
for n = 1:size(indFovVals,2)
    if sum(isnan(indFovVals(:,n)))>0
        indFovVals(:,n) = NaN;
    end
end
    


figure
% plot(steps(2:end)-(mean(diff(steps))/2), indCellVals, 'Color', [0.8 0.8 0.8])
% boxchart(groups(:), indCellVals(:), 'BoxWidth', .03, 'BoxFaceColor', 'k', 'MarkerStyle', '.')
hold on
plotSpread(indFovVals', 'xValues', [0 1], 'distributionColors', 'k')
hold on
plot([0 1], indFovVals, 'k')
errorbar([0 1], nanmean(indFovVals,2), nanstd(indFovVals,[],2)./...
    (sqrt(size(indFovVals,2) - sum(isnan(indFovVals), 2))), 'r')
  ylabel('normalized strength (a.u.)')
    xlabel('V3 actvity')
    ylim([0 2])
% plot(s, movmean(allNrmY(si), .1, 'SamplePoints', s), 'r')

% [N,c] = hist3([allNrmY allNrmX], 'Edges', {0:.01:1, 0:.01:1});

% fittable = table(allNrmX(useCell==1), allNrmY(useCell==1), allCellId(useCell==1), allFishId(useCell==1), 'VariableNames',...
%     {'xVals', 'yVals', 'cellId', 'fishId'});
% lme = fitlme(fittable, 'yVals ~ xVals + (xVals|fishId) + (xVals|fishId:cellId)');
% effectRatio = lme.Coefficients.Estimate(2)/lme.Coefficients.Estimate(1)
[h,p] = ttest(indFovVals(1,:), indFovVals(2,:))
nanmean(indFovVals(2,:)./indFovVals(1,:))

nanmean(indFovVals(1,:))
nanstd(indFovVals(1,:))/sqrt(sum(~isnan(indFovVals(1,:))))
nanmean(indFovVals(2,:))
nanstd(indFovVals(2,:))/sqrt(sum(~isnan(indFovVals(2,:))))
%% Duration vs v3 sr per bout (fig. 4B plot2)
steps = 0:.7:1.5;
allNrmX = [];
allNrmY = [];
finalId = [];
useCell = [];
selectedCells = [data(v3(ismember(v3, find(logSR>0)))).cellId];
selector = ismember(longTable.CellId, selectedCells) & longTable.Freq>freqCutoff(1) &...
    longTable.Freq<freqCutoff(2) & longTable.UseFlag & longTable.SpeedRaw >0.03;
selectedFovs = unique(longTable.FOV(selector));
selectedBouts = unique(longTable.BoutId(selector));

indCellVals = zeros(length(steps)-1,length(selectedCells));
indCellValsSum = zeros(2,length(selectedFovs));

indFovVals = zeros(2, length(selectedFovs));

yVals = longTable.BoutDuration(selector);

xVals = double((longTable.NSpikes(selector)./longTable.BoutDuration(selector))>0);

allCellId = longTable.CellId(selector);
allFishId = longTable.FishId(selector);

for n = 1:length(selectedFovs)
    boutsInFov = longTable.BoutId(selector);
    boutsInFov = boutsInFov(longTable.FOV(selector)==selectedFovs(n));
    

    tmpX = xVals(longTable.FOV(selector)==selectedFovs(n));
    tmpY = yVals(longTable.FOV(selector)==selectedFovs(n));
    

    tmpbin = [];

    nrmX = tmpX;

%     nrmY = tmpY./mean(tmpY);
nrmY = tmpY;
    meanX = zeros(size(unique(boutsInFov)));
    meanY = zeros(size(unique(boutsInFov)));
    bi = 1;
    for k = unique(boutsInFov)'
        meanX(bi) = mean(nrmX(boutsInFov==k));
        meanY(bi) = mean(nrmY(boutsInFov==k));
        bi = bi + 1;
    end
        allNrmX = [allNrmX; meanX];
        allNrmY = [allNrmY; meanY];
        
        for k = 1:2
            indCellValsSum(k,n) = sum(meanX>=(k-1.5) & meanX<(k-0.5));
            if sum(meanX>=k-1.5 & meanX<k-0.5)>1
                indFovVals(k,n) = mean(meanY(meanX>=k-1.5 & meanX<k-0.5));
            else
                indFovVals(k,n) = NaN;

            end
           
        end
end
for n = 1:size(indFovVals,2)
    if sum(isnan(indFovVals(:,n)))>0
        indFovVals(:,n) = NaN;
    end
end
    

figure
% plot(steps(2:end)-(mean(diff(steps))/2), indCellVals, 'Color', [0.8 0.8 0.8])
% boxchart(groups(:), indCellVals(:), 'BoxWidth', .03, 'BoxFaceColor', 'k', 'MarkerStyle', '.')
hold on
plotSpread(indFovVals', 'xValues', [0 1], 'distributionColors', 'k')
hold on
plot([0 1], indFovVals, 'k')
errorbar([0 1], nanmean(indFovVals,2), nanstd(indFovVals,[],2)./...
    (sqrt(size(indFovVals,2) - sum(isnan(indFovVals), 2))), 'r')
  ylabel('duration (s)')
    xlabel('V3 actvity')
    ylim([0 2])

[h,p] = ttest(indFovVals(1,:), indFovVals(2,:))
nanmean(indFovVals(2,:)./indFovVals(1,:))

nanmean(indFovVals(1,:))
nanstd(indFovVals(1,:))/sqrt(sum(~isnan(indFovVals(1,:))))
nanmean(indFovVals(2,:))
nanstd(indFovVals(2,:))/sqrt(sum(~isnan(indFovVals(2,:))))
%% Freq vs v3 SR per bout (fig. 4B plot3)
steps = 0:.7:1.5;
allNrmX = [];
allNrmY = [];
finalId = [];
useCell = [];
selectedCells = [data(v3(ismember(v3, find(logSR>0)))).cellId];
selector = ismember(longTable.CellId, selectedCells) & longTable.Freq>freqCutoff(1) &...
    longTable.Freq<freqCutoff(2) & longTable.UseFlag & longTable.SpeedRaw >0.03;
selectedFovs = unique(longTable.FOV(selector));
selectedBouts = unique(longTable.BoutId(selector));

indCellVals = zeros(length(steps)-1,length(selectedCells));
indCellValsSum = zeros(2,length(selectedFovs));

indFovVals = zeros(2, length(selectedFovs));

yVals = longTable.Freq(selector);

xVals = double((longTable.NSpikes(selector)./longTable.BoutDuration(selector))>0);

allCellId = longTable.CellId(selector);
allFishId = longTable.FishId(selector);
fishId = zeros(length(selectedFovs),1);
for n = 1:length(selectedFovs)
    fishId(n) = unique(allFishId(longTable.FOV(selector)==selectedFovs(n)));
    boutsInFov = longTable.BoutId(selector);
    boutsInFov = boutsInFov(longTable.FOV(selector)==selectedFovs(n));
    

    tmpX = xVals(longTable.FOV(selector)==selectedFovs(n));
    tmpY = yVals(longTable.FOV(selector)==selectedFovs(n));
    

    tmpbin = [];

    nrmX = tmpX;

    nrmY = tmpY;
    meanX = zeros(size(unique(boutsInFov)));
    meanY = zeros(size(unique(boutsInFov)));
    bi = 1;
    for k = unique(boutsInFov)'
        meanX(bi) = mean(nrmX(boutsInFov==k));
        meanY(bi) = mean(nrmY(boutsInFov==k));
        bi = bi + 1;
    end
        allNrmX = [allNrmX; meanX];
        allNrmY = [allNrmY; meanY];
        
        for k = 1:2
            indCellValsSum(k,n) = sum(meanX>=(k-1.5) & meanX<(k-0.5));
            if sum(meanX>=k-1.5 & meanX<k-0.5)>1
                indFovVals(k,n) = mean(meanY(meanX>=k-1.5 & meanX<k-0.5));
            else
                indFovVals(k,n) = NaN;

            end
           
        end
end
for n = 1:size(indFovVals,2)
    if sum(isnan(indFovVals(:,n)))>0
        indFovVals(:,n) = NaN;
    end
end
    

figure
% plot(steps(2:end)-(mean(diff(steps))/2), indCellVals, 'Color', [0.8 0.8 0.8])
% boxchart(groups(:), indCellVals(:), 'BoxWidth', .03, 'BoxFaceColor', 'k', 'MarkerStyle', '.')
hold on
plotSpread(indFovVals', 'xValues', [0 1], 'distributionColors', 'k')
hold on
plot([0 1], indFovVals, 'k')
errorbar([0 1], nanmean(indFovVals,2), nanstd(indFovVals,[],2)./...
    (sqrt(size(indFovVals,2) - sum(isnan(indFovVals), 2))), 'r')
  ylabel('TBF (Hz)')
    xlabel('V3 actvity')
% plot(s, movmean(allNrmY(si), .1, 'SamplePoints', s), 'r')

% [N,c] = hist3([allNrmY allNrmX], 'Edges', {0:.01:1, 0:.01:1});
ylim([0 27])
% fittable = table(allNrmX(useCell==1), allNrmY(useCell==1), allCellId(useCell==1), allFishId(useCell==1), 'VariableNames',...
%     {'xVals', 'yVals', 'cellId', 'fishId'});
% lme = fitlme(fittable, 'yVals ~ xVals + (xVals|fishId) + (xVals|fishId:cellId)');
% effectRatio = lme.Coefficients.Estimate(2)/lme.Coefficients.Estimate(1)
[h,p] = ttest(indFovVals(1,:), indFovVals(2,:))
nanmean(indFovVals(1,:))
nanstd(indFovVals(1,:))/sqrt(sum(~isnan(indFovVals(1,:))))
nanmean(indFovVals(2,:))
nanstd(indFovVals(2,:))/sqrt(sum(~isnan(indFovVals(2,:))))


%% fig4 A examples
figure
tiledlayout(13,1, 'TileSpacing', 'tight', 'Padding', 'tight')


n = find(selectedTrials==98);
idx = 1;
currcell = ismember([data.cellId], selectedCells) & ismember([data.stdV], selectedTrials(n));
currcellidx = find(currcell);
firstCell = currcellidx(1);
yVals = cat(2,data(currcell).trace);
xaxI = xAxImg(500:size(data(find(currcell, 1)).trace,1)+499);
exampleSpikes = double((longTable.NSpikes(selector)./longTable.BoutDuration(selector))>0);
exampleSpikes = exampleSpikes(longTable.recId(selector) == selectedTrials(n));
exampleBoutIds = longTable.BoutId(selector);
exampleBoutIds = exampleBoutIds(longTable.recId(selector) == selectedTrials(n));
exampleCellActive = zeros(size(unique(exampleBoutIds)));
bids = unique(exampleBoutIds);
for k = 1:length(exampleCellActive)
    exampleCellActive(k) = mean(exampleSpikes(exampleBoutIds==bids(k)));
end
srate = length(stdVAll{selectedTrials(n)})/25;
xAxVNR = 1/srate:1/srate:length(stdVAll{selectedTrials(n)})/srate;
exampleBs = data(firstCell).burstVNR(data(firstCell).boutStartIdx);
exampleBe = data(firstCell).burstVNR(data(firstCell).boutEndIdx);
exampleX = mean([exampleBs exampleBe],2);
exampleAmp = longTable.BoutAmp(longTable.CellId == data(firstCell).cellId &...
    longTable.recId == selectedTrials(n));
exampleDist = longTable.SpeedRaw(longTable.CellId == data(firstCell).cellId &...
    longTable.recId == selectedTrials(n));
exampleVNR = stdVAll{selectedTrials(n)};
% exampleOMR = gratingFlow{selectedTrials(n)}.*data(currcellidx(1)).gainAdjust.*...
%     data(currcellidx(1)).omrGain-data(currcellidx(1)).omrSpeed;
exampleOMR = gratingFlow{selectedTrials(n)};
exampleVNR = exampleVNR./mean(exampleAmp);
nexttile(11)
plot(xAxVNR, exampleVNR, 'k')
a1 = gca;
% plot(xAxVNR(exampleBs), (exampleAmp./mean(exampleAmp)).*6-3.4, '.-m')
nexttile(12)
a2 = gca;
hold on
% plot(xAxVNR, exampleOMR.*1.3.*20, 'k');
plot(xAxVNR, exampleOMR, 'k');
hold on
plot(xAxVNR(round(exampleX)), (exampleDist), '.-m')


% yline(-2, 'r')
nexttile(13)
a3 = gca;
plot(xAxVNR(round(exampleX)), (exampleDist), '.-m')
yyaxis right
plot(xAxVNR(round(exampleX)), exampleCellActive, '.-b')

cm = cmocean('-matter', 100);
colorid = exampleAmp./mean(exampleAmp);
colorid = round(nrm(colorid).*99 + 1);
% for bouts = 1:length(exampleAmp)    
%     plot(xAxVNR(exampleBs(bouts):exampleBe(bouts)),...
%         exampleVNR(exampleBs(bouts):exampleBe(bouts)),...
%         'Color', cm(colorid(bouts),:));
% end
nexttile(1,[10 1])
a4 = gca;
hold on
cm = cmocean('-matter', 35);
for k = 2:size(yVals,2)
    if ~isnan(yVals(1,k))
        xs = data(firstCell).burstImg(data(firstCell).boutStartIdx);
        xe = data(firstCell).burstImg(data(firstCell).boutEndIdx);
        
        plot(xaxI, yVals(:,k)-medfilt1(yVals(:,k), 201)+idx*4.5 +3, 'Color', [.3 .3 .3]);
        
        exampleNS = longTable.NSpikes(longTable.CellId == data(currcellidx(k)).cellId &...
            longTable.recId == selectedTrials(n));
        exampleBD = longTable.BoutDuration(longTable.CellId == data(currcellidx(k)).cellId &...
            longTable.recId == selectedTrials(n));
        exampleSR = exampleNS./exampleBD;
%         plot(xaxI(xs-500), (exampleSR./mean(exampleSR))+idx*4.5+3, '.-m')
%         for bouts = 1:length(exampleSR)
%             tmpXidx = (xs(bouts):xe(bouts))-500;
%             tmpXidx(tmpXidx<1) = [];
%             plot(xaxI(tmpXidx), yVals(tmpXidx,k)+idx*5 +3,...
%                 'Color', cm(round(exampleSR(bouts))+1,:));
%         end
        idx = idx + 1;
    end
end
linkaxes([a1 a2 a3 a4], 'x')
xlim([4 10.2])
  
%% spike triggered correlations maps example
%carefull this is a 7Gb file
mov = nrrdread('.\helper_data\rec00001.nrrd');
mov = mov(:,:,501:end);
avFrame = loadtiff('.\helper_data\rec00001_meanFrame.tif');
%%
selector = false(size(data));
for n = 1:length(data)
    selector(n) = data(n).trial==1 &&...
        data(n).fishId==4 &&...
        data(n).zPos==1 &&...
        ~isempty(data(n).av);
end

cId1 = cat(1, data(selector).cellId);

selector = false(size(data));
for n = 1:length(data)
    selector(n) = data(n).trial==1 &&...
        data(n).fishId==4 &&...
        data(n).zPos==1 &&...
        ~isempty(data(n).av);
end

cId2 = cat(1, data(selector).cellId);

cId = ismember(cId1, cId2);

cId = cId1(cId);

selector = false(size(data));
for n = 1:length(data)
    selector(n) = (data(n).trial==1 ||...
        data(n).trial==1) &&...
        ismember(data(n).cellId, cId);
end


av = cat(2,data(selector).av);
amp = max(av) - min(av);
ids = find(selector);
% ids = ids(amp>=.5);

tr = cat(2,data(ids).trace);
xPos = cat(1,data(ids).xyzpos);
xPos = xPos(:,1);
yPos = cat(1,data(ids).xyzpos);
yPos = yPos(:,2);
spikes = cat(1,{data(ids).spikes});
% Figure S5, D,E
spikeTimeMat = zeros(length(spikes));
spikeTimeMat2 = zeros(length(spikes));
for d1 = 1:length(spikes)
    for d2 = 1:length(spikes)
        tmpSpikes = 0;
        for nSpikes = 1:length(spikes{d1})
            tmpSpikes = tmpSpikes + sum(ismember(spikes{d1}(nSpikes)-2:spikes{d1}(nSpikes)+2,...
                spikes{d2}));
        end
        spikeTimeMat(d1, d2) = tmpSpikes;
    end
    spikeTimeMat2(d1,:) = spikeTimeMat(d1,:)./spikeTimeMat(d1,d1);
end
figure
h = heatmap(round(spikeTimeMat,2));
cmocean('thermal')
h.GridVisible = 'off';
figure
h = heatmap(round(spikeTimeMat2,2));
cmocean('thermal')
h.GridVisible = 'off';

windowSize = 40;
corrImgStack = zeros(size(mov,1), size(mov,2), length(spikes));
for k = 1:length(spikes)
stMov = zeros(size(mov,1), size(mov,2), windowSize + 1);
cellN = k;
tmpSp = spikes{cellN};
stTrace = zeros(windowSize + 1,length(tmpSp));

% get first spike estimate
for n = 1:length(tmpSp)
    if (tmpSp(n)-windowSize/2)>0 && tmpSp(n)+windowSize/2<=length(tr(:,cellN))
        stTrace(:,n) = tr(tmpSp(n)-windowSize/2:tmpSp(n)+windowSize/2,cellN);
    end
end
meanSpike = mean(stTrace-mean(stTrace),2);

spxc = xcorr(tr(:,cellN), meanSpike);
spxc = spxc(size(tr,1)-windowSize/2:end);

%refine spike timing with xCorr
pind = zeros(size(tmpSp));
for n = 1:length(tmpSp)
    [~, pind(n)] = max(spxc(tmpSp(n)-2:tmpSp(n)+2));
    pind(n) = pind(n) - 3;
end
pind = pind+tmpSp;

for n = 1:length(tmpSp)
    if (tmpSp(n)-windowSize/2)>0 && tmpSp(n)+windowSize/2<=length(tr(:,cellN)) &&...
            spxc(pind(n))>(mean(spxc(pind)) - 2*std(spxc(pind)))
%         tmpMov = double(mov(:,:,tmpSp(n)-windowSize/2:tmpSp(n)+windowSize/2));
        tmpMov = double(mov(:,:,pind(n)-windowSize/2:pind(n)+windowSize/2));
        tmpMov = tmpMov-mean(tmpMov(:));
        stMov = stMov + tmpMov./length(tmpSp);
        stTrace(:,n) = tr(pind(n)-windowSize/2:pind(n)+windowSize/2,cellN);
%         stTrace(:,n) = tr(tmpSp(n)-windowSize/2:tmpSp(n)+windowSize/2,cellN);
    end
end
meanSpike = mean(stTrace-mean(stTrace),2);
dStMov = stMov - repmat(mean(stMov,3), [1 1 size(stMov,3)]);
% saveastiff(single(dStMov), ['cell_' num2str(k) '.tif']);
% for n = 1:size(dStMov,3)
%     dStMov(:,:,n) = imgaussfilt(dStMov(:,:,n),1);
%     stMov(:,:,n) = imgaussfilt(stMov(:,:,n),1);
% end
corrImg = corr(tovec(dStMov)', meanSpike);
corrImgStack(:,:,k) = reshape(corrImg, size(stMov,1), size(stMov,2));
end
%%
corrImgStackfilt = zeros(size(corrImgStack));
for n = 1:size(corrImgStack, 3)
    corrImgStackfilt(:,:,n) = medfilt2(corrImgStack(:,:,n), [4 4]);
end
%% figure 3E
figure
tiledlayout(size(corrImgStackfilt,3)/2, 2, 'Padding', 'none', 'Tilespacing', 'none')
cm = distinguishable_colors(size(corrImgStackfilt,3), [0 0 0]);
fullImg = zeros(98, 932, 3);
% montage = zeros(98*11, 932*2);
for n = 1:size(corrImgStackfilt, 3)
    nexttile
    tI = corrImgStackfilt(:,:,n);
    tI = imrotate(tI, 2, 'bicubic', 'crop');
    tI = tI(23:23+97,39:39+931);
    tI(tI<0) = 0;
    tI(tI>.95) = .95;
    tI = tI./0.95;
%     montage(1:size(tI,1), 1:size(tI,2)) = tI;
    tmpImg = zeros(size(tI,1), size(tI,2), 3);
    imagesc(tI);
%     saveastiff(tI, ['X:\Lab\Papers\zf_spinal_imaging\Fig6\cell' num2str(n) '.tif']);
    colormap('gray')
%     tI = 1-tI;
%     r = ceil(rand.*100);
    tmpImg(:,:,1) = tI.*cm(n,1);
    tmpImg(:,:,2) = tI.*cm(n,2);
    tmpImg(:,:,3) = tI.*cm(n,3);
    fullImg = max(cat(4, fullImg, tmpImg),[],4);
%     colormap('gray');
%     imshow(tmpImg)
    daspect([1 1 1])
    axis off
end
% saveastiff(fullImg, 'X:\Lab\Papers\zf_spinal_imaging\Fig6\colored_composit.tif');
figure
imshow(fullImg)

avCrop = imrotate(avFrame, 2, 'bicubic', 'crop');
avCrop = avCrop(23:23+97,39:39+931);
% saveastiff(avCrop, 'X:\Lab\Papers\zf_spinal_imaging\Fig6\avFrame_croped.tif')
%% SR average per cell and suthreshold average per cell

% use this for Fig 2A and supl fig
sel = v3(ismember(v3, find(logSR>0)));

% don't use this...use this for Supl. Fig. 1C and 1D
% sel = nonOsc;

trRemS = zeros(size(tr,1), length(sel));
for n = 1:length(sel)
    tmpS = data(sel(n)).spikes;
    tmpTr = data(sel(n)).trace;
    for k = 1:length(tmpS)
        p1 = max([tmpS(k)-3 1]);
        p2 = min([tmpS(k)+3 length(tmpTr )]);
        for p = p1:p2
            tmpTr(p) = median(tmpTr(max([1 p-5]):p));
        end
    end
    trRemS(:,n) = tmpTr;
end
tr = [data(sel).trace];

subthrV = [];
RawV = [];
RawVend = [];
subthrVend = [];
srStart = [];
srEnd = [];
speedStart = [];
speedEnd = [];
vnrStart = [];
vnrEnd = [];
cellId = [];
idx = 1;
for n = 1:length(sel)
    disp(num2str(n))
    ba = data(sel(n)).burstImg-499;
    bV = data(sel(n)).burstVNR;
        bs = ba(data(sel(n)).boutStartIdx);
        be = ba(data(sel(n)).boutEndIdx);
        bsV = bV(data(sel(n)).boutStartIdx);
        beV = bV(data(sel(n)).boutEndIdx);
        %get rid of bouts that are outside the imaging time
        noBout = bs<1 | be<1 | bs>size(tr,1) | be>size(tr,1);
        bs(noBout) = [];
        be(noBout) = [];
        bsV(noBout) = [];
        beV(noBout) = [];
        noBurst = ba<bs(1) | ba>be(end);
        ba(noBurst) = [];
        currSpeed = gratingFlow{data(sel(n)).stdV};
        currSpeed = currSpeed.*data(sel(n)).gainAdjust;
        currVNR = stdVAll{data(sel(n)).stdV};
        vfact = length(currSpeed)/25/1000;
    for k  = 1:length(bs)
        tmpTr = nan(401,1);
        tmpTrRaw = nan(401,1);
        tmpSR = zeros(401,1);
        tmpSpeed = zeros(401*vfact,1);
        tmpVNR = zeros(401*vfact,1);
        %always deal with edge cases
        p1 = min([200 bs(k)-1]);
        p2 = min([200 size(trRemS,1)-bs(k)]);
       
        tmpTr(201-p1:201+p2) = trRemS(bs(k)-p1:bs(k)+p2,n);
        tmpTrRaw(201-p1:201+p2) = tr(bs(k)-p1:bs(k)+p2,n);
        triggeredSpikes = data(sel(n)).spikes - bs(k) + 201;
        triggeredSpikes = triggeredSpikes(triggeredSpikes>0 & triggeredSpikes<=401);
        tmpSR(triggeredSpikes) = 1;
        p1V = min([200*vfact bsV(k)-1]);
        p2V = min([200*vfact length(currSpeed)-bsV(k)]);
        tmpSpeed(201*vfact-p1V:201*vfact+p2V) = currSpeed(bsV(k)-p1V:bsV(k)+p2V);
        tmpVNR(201*vfact-p1V:201*vfact+p2V) = currVNR(bsV(k)-p1V:bsV(k)+p2V);
        %         subthrV = cat(2, subthrV, tmpTr-mean(tmpTr, 'omitnan'));
        F0 = tmpTr(find(~isnan(tmpTr),1));
        subthrV = cat(2, subthrV, tmpTr-F0);
        RawV = cat(2, RawV, tmpTrRaw);
        srStart = cat(2, srStart, tmpSR);
        if length(tmpSpeed)==4010
            tmpSpeed = interp1(1:4010, tmpSpeed, .2:.2:4010);
            tmpSpeed = tmpSpeed';
            tmpVNR = interp1(1:4010, tmpVNR, .2:.2:4010);
            tmpVNR = tmpVNR';
        end
        speedStart = cat(2, speedStart, tmpSpeed);
        vnrStart = cat(2, vnrStart, tmpVNR);
        % subthrV = cat(2, subthrV, nrm(tmpTr));
        
        tmpTr = nan(401,1);
        tmpTrRaw = nan(401,1);
        tmpSR = zeros(401,1);
        tmpSpeed = zeros(401*vfact,1);
        tmpVNR = zeros(401*vfact,1);
        p1 = min([200 be(k)-1]);
        p2 = min([200 size(trRemS,1)-be(k)]);
        tmpTr(201-p1:201+p2) = trRemS(be(k)-p1:be(k)+p2,n);
        tmpTrRaw(201-p1:201+p2) = tr(be(k)-p1:be(k)+p2,n);
        triggeredSpikes = data(sel(n)).spikes - be(k) + 201;
        triggeredSpikes = triggeredSpikes(triggeredSpikes>0 & triggeredSpikes<=401);
        tmpSR(triggeredSpikes) = 1;
        p1V = min([200*vfact beV(k)-1]);
        p2V = min([200*vfact length(currSpeed)-beV(k)]);
        tmpSpeed(201*vfact-p1V:201*vfact+p2V) = currSpeed(beV(k)-p1V:beV(k)+p2V);
        tmpVNR(201*vfact-p1V:201*vfact+p2V) = currVNR(beV(k)-p1V:beV(k)+p2V);
        %         subthrVend = cat(2, subthrVend, tmpTr-mean(tmpTr, 'omitnan'));
        subthrVend = cat(2, subthrVend, tmpTr-F0);
        RawVend = cat(2, RawVend, tmpTrRaw);
        srEnd = cat(2, srEnd, tmpSR);
        if length(tmpSpeed)==4010
            tmpSpeed = interp1(1:4010, tmpSpeed, .2:.2:4010);
            tmpSpeed = tmpSpeed';
            tmpVNR = interp1(1:4010, tmpVNR, .2:.2:4010);
            tmpVNR = tmpVNR';
        end
        speedEnd = cat(2, speedEnd, tmpSpeed);
        vnrEnd = cat(2, vnrEnd, tmpVNR);
        % subthrVend = cat(2, subthrVend, nrm(tmpTr));
        cellId(idx) = sel(n);
        idx = idx + 1;
    end
end
% [~, si] = sort(logSR(cellId));
% imagesc(subthrV(:,si)')
avsubthV = zeros(401,length(sel));
avsubthVend = zeros(401,length(sel));
avSrStart = zeros(401,length(sel));
avSrEnd = zeros(401,length(sel));

for n = 1:size(avsubthV, 2)
    avsubthV(:,n) = mean(subthrV(:,ismember(cellId, sel(n))),2, 'omitnan');
    avsubthVend(:,n) = mean(subthrVend(:,ismember(cellId, sel(n))),2, 'omitnan');
    avSrStart(:,n) = mean(srStart(:,ismember(cellId, sel(n))),2);
    avSrStart(:,n) = movmean(avSrStart(:,n), 20);
    avSrEnd(:,n) = mean(srEnd(:,ismember(cellId, sel(n))),2);
    avSrEnd(:,n) = movmean(avSrEnd(:,n), 20);
end
trials = unique([data(sel).stdV]);
avSpStart = zeros(size(speedStart,1),length(trials));
avSpEnd = zeros(size(speedEnd,1),length(trials));
avVnrStart = zeros(size(speedStart,1),length(trials));
avVnrEnd = zeros(size(speedEnd,1),length(trials));
nbtstot = 0;
for n = 1:length(trials) 
    avSpStart(:,n) = mean(speedStart(:,ismember([data(cellId).stdV], trials(n))),2, 'omitnan');
    avSpEnd(:,n) = mean(speedEnd(:,ismember([data(cellId).stdV], trials(n))),2, 'omitnan');
    avVnrStart(:,n) = mean(vnrStart(:,ismember([data(cellId).stdV], trials(n))),2, 'omitnan');
    avVnrEnd(:,n) = mean(vnrEnd(:,ismember([data(cellId).stdV], trials(n))),2, 'omitnan');
    nbtstot = nbtstot + size([data([data.stdV] == trials(n)).boutStartIdx],1);
end
nbts = zeros(size(sel));
for n = 1:length(sel)
    nbts(n) = sum(cellId==sel(n));
end

%% Fig. 2A
xa = -200:200;
figure
tiledlayout(5,2,'Padding', 'none')
nexttile(1, [2 1])
[~, si] = sort(logSR(sel), 'descend');
imagesc([-200 200], [1 size(avSrStart,2)], avSrStart(:,si)', [min([avSrStart(:); avSrEnd(:)]) max([avSrStart(:); avSrEnd(:)])])
cmocean('-matter')
xline(0,'b--', 'LineWidth', 2)
xlim([-100 200])
ylabel('cell #')
nexttile(5)
av = mean(avSrStart(:,logSR(sel)>0),2, 'omitnan');
sem = std(avSrStart(:,logSR(sel)>0), [], 2, 'omitnan')./sqrt(size(avSrStart(:,logSR(sel)>0),2));
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 -Inf];
% [fob, gof] = fit(xa(find(xa==-37):find(xa==149))', av(find(xa==-37):find(xa==149)), 'c-a*exp(-x/b)', opts);
[fob, gof] = fit(xa(find(xa==-37):find(xa==40))', av(find(xa==-37):find(xa==40)), 'a*x+b', opts);
tau_start = fob.a;
alpha = 0.95;
ci = confint(fob, alpha); 
t = tinv((1+alpha)/2, gof.dfe); 
se = (ci(2,2)-ci(1,2)) ./ (2*t); % Standard Error
disp(['\tau start = ' num2str(tau_start) ' ± ' num2str(se) ' ms'])
plot(-200:200, av, 'k')
patch([-200:200 200:-1:-200], [av+sem; flipud(av-sem)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on
plot(-37:40,fob([-37:40]),'r')
ylim([0 0.011])
xlim([-100 200])
ylabel('frequency (kHz)')
nexttile(7)
av = mean(avSrStart(:,logSR(sel)<=0),2, 'omitnan');
sem = std(avSrStart(:,logSR(sel)<=0), [], 2, 'omitnan')./sqrt(size(avSrStart(:,logSR(sel)<=0),2));
plot(-200:200, av, 'k')
patch([-200:200 200:-1:-200], [av+sem; flipud(av-sem)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
ylim([0 0.011])
xlim([-100 200])
nexttile(9)
hold on
ylabel('VNR')
av = mean(avSpStart,2, 'omitnan');
sem = std(avSpStart, [], 2, 'omitnan')./sqrt(size(avSpStart,2));
plot((1/50:1/50:401)-201, av, 'k')
patch([(1/50:1/50:401)-201 (401:-1/50:1/50)-201], [av+sem; flipud(av-sem)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

av = mean(avVnrStart,2, 'omitnan')-0.005;
sem = std(avVnrStart, [], 2, 'omitnan')./sqrt(size(avVnrEnd,2));
plot((1/50:1/50:401)-201, av.*150, 'b')
patch([(1/50:1/50:401)-201 (401:-1/50:1/50)-201], [av+sem; flipud(av-sem)].*150, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
ylim([-0.1 1.2])
xlim([-100 200])

nexttile(2, [2 1])
imagesc([-200 200], [1 size(avSrEnd,2)], avSrEnd(:,si)', [min([avSrStart(:); avSrEnd(:)]) max([avSrStart(:); avSrEnd(:)])])
cmocean('-matter')
xline(0,'b--', 'LineWidth', 2)
xlim([-100 200])
nexttile(6)
av = mean(avSrEnd(:,logSR(sel)>0),2, 'omitnan');
sem = std(avSrEnd(:,logSR(sel)>0), [], 2, 'omitnan')./sqrt(size(avSrEnd(:,logSR(sel)>0),2));
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf];
[fob, gof] = fit(xa(find(xa==0):find(xa==183))', av(find(xa==0):find(xa==183)), 'a*x+b', opts);
tau_end = fob.a;
alpha = 0.95;
ci = confint(fob, alpha); 
t = tinv((1+alpha)/2, gof.dfe); 
se = (ci(2,2)-ci(1,2)) ./ (2*t); % Standard Error
disp(['\tau end = ' num2str(tau_end) ' ± ' num2str(se) ' ms'])

plot(-200:200, av, 'k')
patch([-200:200 200:-1:-200], [av+sem; flipud(av-sem)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on
plot(0:183,fob([0:183]),'r')
ylim([0 0.011])
xlim([-100 200])
nexttile(8)
av = mean(avSrEnd(:,logSR(sel)<=0),2, 'omitnan');
sem = std(avSrEnd(:,logSR(sel)<=0), [], 2, 'omitnan')./sqrt(size(avSrEnd(:,logSR(sel)<=0),2));
plot(-200:200, av, 'k')
patch([-200:200 200:-1:-200], [av+sem; flipud(av-sem)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
ylim([0 0.011])
xlim([-100 200])

nexttile(10)
hold on

av = mean(avSpEnd,2, 'omitnan');
sem = std(avSpEnd, [], 2, 'omitnan')./sqrt(size(avSpEnd,2));
plot((1/50:1/50:401)-201, av, 'k')
patch([(1/50:1/50:401)-201 (401:-1/50:1/50)-201], [av+sem; flipud(av-sem)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
av = mean(avVnrEnd,2, 'omitnan')-0.005;
sem = std(avVnrEnd, [], 2, 'omitnan')./sqrt(size(avVnrEnd,2));
plot((1/50:1/50:401)-201, av.*150, 'b')
patch([(1/50:1/50:401)-201 (401:-1/50:1/50)-201], [av+sem; flipud(av-sem)].*150, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
ylim([-0.1 1.2])
xlim([-100 200])


nexttile(7)
av = mean(avsubthV(:,logSR(sel)>0),2, 'omitnan');
sem = std(avsubthV(:,logSR(sel)>0), [], 2, 'omitnan')./sqrt(size(avsubthV(:,logSR(sel)>0),2));
plot(-200:200, av, 'k')
patch([-200:200 200:-1:-200], [av+sem; flipud(av-sem)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
ylim([-0.2 0.8])
xlim([-100 200])
ylabel('subthr. fluo.')

nexttile(8)
av = mean(avsubthVend(:,logSR(sel)>0),2, 'omitnan');
sem = std(avsubthVend(:,logSR(sel)>0), [], 2, 'omitnan')./sqrt(size(avsubthVend(:,logSR(sel)>0),2));
plot(-200:200, av, 'k')
patch([-200:200 200:-1:-200], [av+sem; flipud(av-sem)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
ylim([-.2 0.8])
xlim([-100 200])

figure
bar = linspace(min([avSrStart(:); avSrEnd(:)]), max([avSrStart(:); avSrEnd(:)]),100);
bar = repmat(bar,10,1);
imagesc([min([avSrStart(:); avSrEnd(:)]) max([avSrStart(:); avSrEnd(:)])].*1000, [0 1], bar)
cmocean('-matter')

%% number of spikes/recording in nonOsc cells
nSpikes = zeros(size(nonOsc));
for n  = 1:length(nonOsc)
    nSpikes(n) = length(data(nonOsc(n)).spikes);
end
%% suppl mov 1
selector = false(size(data));
for n = 1:length(data)
    selector(n) = data(n).trial==1 &&...
        data(n).fishId==4 &&...
        data(n).zPos==1 &&...
        ~isempty(data(n).av);
end

cId1 = cat(1, data(selector).cellId);

selector = false(size(data));
for n = 1:length(data)
    selector(n) = data(n).trial==1 &&...
        data(n).fishId==4 &&...
        data(n).zPos==1 &&...
        ~isempty(data(n).av);
end

cId2 = cat(1, data(selector).cellId);

cId = ismember(cId1, cId2);

cId = cId1(cId);

selector = false(size(data));
for n = 1:length(data)
    selector(n) = (data(n).trial==1 ||...
        data(n).trial==1) &&...
        ismember(data(n).cellId, cId);
end


av = cat(2,data(selector).av);
amp = max(av) - min(av);
ids = find(selector);
% ids = ids(amp>=.5);

tr = cat(2,data(ids).trace);
stdv = stdVAll{data(ids(1)).stdV};
speed = gratingFlow{data(ids(1)).stdV};

srImg = 996.3;
xAxImg = 1/srImg:1/srImg:25000/srImg;
srate = 50000;
xAxVNR = 1/srate:1/srate:length(stdVAll{1})/srate;
% clearvars -except xAxImg xAxVNR stdv speed tr
%% supplementary movie 1
%carefull this is a 7Gb file
% mov = nrrdread('C:\Users\Urs\Desktop\rec00001.nrrd');
mov = nrrdread('X:\Lab\Labmembers\Urs Boehm\Data\7-zArchon\20190228\l2\rec00001.nrrd');

% avFrame = loadtiff('X:\Lab\Labmembers\Urs Boehm\Data\7-zArchon\20190228\l2\rec00001_meanFrame.tif');
% mov = loadtiff("C:\Users\Urs\Documents\substack3to5s.tif");
%% define time snipet and calculate df movie
frames_start = find(xAxImg>3.2,1);
frames_end = find(xAxImg<4.5,1,'last');
mov = mov(:,:,frames_start:frames_end);
avFrame = mean(mov,3);
dfmov = double(mov)-avFrame;
img3 = 0.8*repmat(mat2gray(avFrame), [1 1 3]);


%% rotate and crop movie and denoise with PCA
dfr = imrotate3(dfmov,2.0,[0 0 1]);
dfrc = dfr(52:127,15:1015,:);
[ysize, xsize, ~] = size(dfrc);
[a, b, c] = pca(tovec(dfrc));
pcaKeep = 30;
res = reshape(b(:,1:pcaKeep)*a(:,1:pcaKeep)',ysize, xsize, []);
%% generate color movie
avr = imrotate(avFrame,2.0);
avrc = avr(52:127,15:1015);
img3 = 0.8*repmat(mat2gray(avrc), [1 1 3]);
cmin4 = 10;
cmax4 = 80;
Colormov = zeros(size(res,1), size(res,2), 3, size(res,3));
figure
clear M
for j = 1:size(Colormov,4)
    % .8*img3 + 
    Colormov(:,:,:,j) = .8*img3 + .8*grs2rgb(res(:,:,j), colormap(hot), 5.8, 20.2);
    imshow(Colormov(:,:,:,j))
    text(135,15,[num2str(xAxImg(frames_start + j), '%+3.3f') ' s'], 'FontSize', 10, 'color', [1 1 1], 'HorizontalAlignment', 'right')
    drawnow
    M(j) = getframe(gca);
end
%% save pca denoised color movie
%     v = VideoWriter(['C:\Users\Urs\Documents\pca_color.avi'], 'Uncompressed AVI');
%     open(v)
%     writeVideo(v,M)
%     close(v)
%% color movie with traces underneath
    nToDisp = 5;
    roid = find(selector);
    roid = roid([1 7 11 14 20]);
    rois = [data(roid).rois];
    centerA = fliplr(size(avFrame)./2);
    centerB = fliplr(size(avr)./2);

    for k = 1:length(rois)
        %rotate
        rois{k} = (rois{k}-centerA)*[cosd(2.0) -sind(2.0); sind(2.0) cosd(2.0)]+centerB;
        %crop
        rois{k}(:,2) = rois{k}(:,2)-51;
        rois{k}(:,1) = rois{k}(:,1)-14;
    end
    f = figure;
    f.Position = [50 50 942 453];
    set(f, 'Color', [0 0 0])
    
    tiledlayout(6,1, 'Padding', 'tight', 'TileSpacing', 'none')
    
    % nexttile
    % imshow2(avFrame,[])
    % nexttile([2 1])
    xax = xAxImg(501:501+length(tr)-1);
    nexttile(1)
    imshow(M(1).cdata)
   
    nexttile([4 1])
    stackplot(tr(:,1:nToDisp),xax)
    xlim([3.2 4.5])
    nexttile(6)
    plot(xAxVNR, stdv, 'Color', [.5 .5 .5])
    
    xlim([3.2 4.5])
 clear MM   
    for j = 1:(length(xax(frames_start:frames_end)))
        nexttile(1)
        imshow(Colormov(:,:,:,j))
        hold on
        text(110,10,[num2str(xAxImg(frames_start + j), '%+3.3f') ' s'], 'FontSize', 20, 'color', [1 1 1], 'HorizontalAlignment', 'right')
        for k = 1:length(rois)
            plot(rois{k}(:,1), rois{k}(:,2))
        end
        hold off
        nexttile(2)
        stackplot(tr(:,[1 7 11 14 20]),xax)
        set(gca,'Color', [0 0 0])
        hold on
        xline(xax(find(xax>3.2,1) + j), 'w')
        xlim([3.2 4.5])
        xlabel('time (s)')
        ylabel('fluorescence (a.u.)')
        text(3.18,4, 'fluorescence', 'Rotation', 90, 'Color', 'w')
        set(gca,'visible', 'off')
        
%         set(gca,'xcolor','w')
%         set(gca,'ycolor','w') 
%         set(gca,'YTickLabel',[])
%         set(gca,'YTick',[])
        hold off
        nexttile(6)
        set(gca,'Color', [0 0 0])
        plot(xAxVNR, stdv, 'Color', [.5 .5 .5])
        set(gca,'Color', [0 0 0])
        hold on
        xline(xax(find(xax>3.2,1) + j), 'w')
        xlim([3.2 4.5])
        xlabel('time (s)')
        ylabel('VNR')
        set(gca,'box', 'off')
        set(gca,'xcolor','w')
        set(gca,'ycolor','w') 
        set(gca,'YTickLabel',[])
        set(gca,'YTick',[])
         hold off
        drawnow
       
        MM(j) = getframe(gcf);
    end

%% Fig. S2 plot a few high snr V2as
cids = [1, 5, 27, 34, 36, 191, 204];
% figure
% stackplot(tr(:,cids), xAxImg(501:501+length(tr)-1))
figure
% tiledlayout(length(unique([data(V2as(cids)).stdV])), 1,...
%     'TileSpacing', 'tight', 'Padding', 'tight')
currstdv = 0;
xa = xAxImg(501:501+length(tr)-1);
shiftn = 0;
hold on
for n = 1:length(cids)
    if data(V2as(cids(n))).stdV ~= currstdv
        currstdv = data(V2as(cids(n))).stdV;
        srate = length(stdVAll{data(V2as(cids(n))).stdV})/25;
        xAxVNR = 1/srate:1/srate:length(stdVAll{data(V2as(cids(n))).stdV})/srate;
        plot(xAxVNR(xAxVNR>.5 & xAxVNR < 3.5),...
            nrm(stdVAll{data(V2as(cids(n))).stdV}(xAxVNR>.5 & xAxVNR < 3.5)) + shiftn*1.1, 'k')
        shiftn = shiftn + 1;
    end
    plot(xa(xa>0.5 & xa<3.5),...
        nrm(data(V2as(cids(n))).trace(xa>0.5 & xa<3.5))+ shiftn*1.1)
    
    shiftn = shiftn + 1;
end
%% suppl mov 2
selector = false(size(data));
for n = 1:length(data)
    selector(n) = data(n).trial==1 &&...
        data(n).fishId==4 &&...
        data(n).zPos==2 &&...
        ~isempty(data(n).av);
end

cId1 = cat(1, data(selector).cellId);

selector = false(size(data));
for n = 1:length(data)
    selector(n) = data(n).trial==1 &&...
        data(n).fishId==4 &&...
        data(n).zPos==2 &&...
        ~isempty(data(n).av);
end

cId2 = cat(1, data(selector).cellId);

cId = ismember(cId1, cId2);

cId = cId1(cId);

selector = false(size(data));
for n = 1:length(data)
    selector(n) = (data(n).trial==1 ||...
        data(n).trial==1) &&...
        ismember(data(n).cellId, cId);
end


av = cat(2,data(selector).av);
amp = max(av) - min(av);
ids = find(selector);
% ids = ids(amp>=.5);

tr = cat(2,data(ids).trace);
stdv = stdVAll{data(ids(1)).stdV};
speed = gratingFlow{data(ids(1)).stdV};

srImg = 996.3;
xAxImg = 1/srImg:1/srImg:25000/srImg;
srate = 50000;
xAxVNR = 1/srate:1/srate:length(stdVAll{1})/srate;
% clearvars -except xAxImg xAxVNR stdv speed tr
%% supplementary movie 2
%carefull this is a 7Gb file
mov = nrrdread(".\helper_data\rec00005.nrrd");
% mov = nrrdread("C:\Users\Urs\Desktop\rec00005.nrrd");
% mov = nrrdread('X:\Lab\Labmembers\Urs Boehm\Data\7-zArchon\20190228\l2\rec00001.nrrd');

% avFrame = loadtiff('X:\Lab\Labmembers\Urs Boehm\Data\7-zArchon\20190228\l2\rec00001_meanFrame.tif');
% mov = loadtiff("C:\Users\Urs\Documents\substack3to5s.tif");
%% define time snipet and calculate df movie
frames_start = find(xAxImg>3,1);
frames_end = find(xAxImg<4,1,'last');
mov = mov(:,:,frames_start:frames_end);
avFrame = mean(mov,3);
dfmov = double(mov)-avFrame;
img3 = 0.8*repmat(mat2gray(avFrame), [1 1 3]);

% %% put time stamp and save df movie
% for j = 1:size(Colormov,4)
%     imshow(dfmov(:,:,j), [-100 100])
%     text(110,10,[num2str(xAxImg(frames_start + j), '%+3.3f') ' s'], 'FontSize', 20, 'color', [1 1 1], 'HorizontalAlignment', 'right')
%     drawnow
%     M(j) = getframe(gca);
% end
% v = VideoWriter(['C:\Users\Urs\Desktop\df_movie_gray.avi'], 'Uncompressed AVI');
%     open(v)
%     writeVideo(v,M)
%     close(v)
%     
%% rotate and crop movie and denoise with PCA
dfr = imrotate3(dfmov,2.0,[0 0 1]);
dfrc = dfr(40:137,:,:);
[ysize, xsize, ~] = size(dfrc);
[a, b, c, ~, explained] = pca(tovec(dfrc));
pcaKeep = 6;
res = reshape(b(:,1:pcaKeep)*a(:,1:pcaKeep)',ysize, xsize, []);
%% generate color movie
avr = imrotate(avFrame,2.0);
avrc = avr(40:137,:);
img3 = 0.8*repmat(mat2gray(avrc), [1 1 3]);
cmin4 = 10;
cmax4 = 80;
Colormov = zeros(size(res,1), size(res,2), 3, size(res,3));
figure
clear M
for j = 1:size(Colormov,4)
    % .8*img3 + 
    Colormov(:,:,:,j) = .8*img3 + .8*grs2rgb(res(:,:,j), colormap(hot), 3, 10);
    imshow(Colormov(:,:,:,j))
    text(135,15,[num2str(xAxImg(frames_start + j), '%+3.3f') ' s'], 'FontSize', 10, 'color', [1 1 1], 'HorizontalAlignment', 'right')
    drawnow
    M(j) = getframe(gca);
end
% %% save pca denoised color movie
%     v = VideoWriter(['X:\Lab\Labmembers\Urs Boehm\Data\7-zArchon\pca_color.avi'], 'Uncompressed AVI');
%     open(v)
%     writeVideo(v,M)
%     close(v)
%% color movie with traces underneath
    nToDisp = 4;
    roid = find(selector);
    roid = roid([7 8 18 23]);
    rois = [data(roid).rois];
    centerA = fliplr(size(avFrame)./2);
    centerB = fliplr(size(avr)./2);

    for k = 1:length(rois)
        %rotate
        rois{k} = (rois{k}-centerA)*[cosd(2.0) -sind(2.0); sind(2.0) cosd(2.0)]+centerB;
        %crop
        rois{k}(:,2) = rois{k}(:,2)-40;
        rois{k}(:,1) = rois{k}(:,1);
    end
    f = figure;
    f.Position = [145.8000   84.6000  589.2000  484.0000];
    set(f, 'Color', [0 0 0])
    
    tiledlayout(6,1, 'Padding', 'compact', 'TileSpacing', 'none')
    
    % nexttile
    % imshow2(avFrame,[])
    % nexttile([2 1])
    xax = xAxImg(501:501+length(tr)-1);
    nexttile(1)
    imshow(M(1).cdata)
   
    nexttile([4 1])
    stackplot(tr(:,1:nToDisp),xax)
    xlim([3 4])
    nexttile(6)
    plot(xAxVNR, stdv, 'Color', [.5 .5 .5])
    
    xlim([3 4])
 clear MM   
    for j = 1:(length(xax(frames_start:frames_end)))
        nexttile(1)
        imshow(Colormov(:,:,:,j))
        hold on
        text(135,15,[num2str(xAxImg(frames_start + j), '%+3.3f') ' s'], 'FontSize', 10, 'color', [1 1 1], 'HorizontalAlignment', 'right')
        for k = 1:length(rois)
            plot(rois{k}(:,1), rois{k}(:,2))
        end
        hold off
        nexttile(2)
        stackplot(tr(:,[7 8 18 23]),xax)
        set(gca,'Color', [0 0 0])
        hold on
        xline(xax(find(xax>3,1) + j), 'w')
        xlim([3 4])
        xlabel('time (s)')
        ylabel('fluorescence (a.u.)')
        
        text(2.981071155317521,2.740647482014388, 'fluorescence', 'Rotation', 90, 'Color', 'w')
        set(gca,'visible', 'off')
        
%         set(gca,'xcolor','w')
%         set(gca,'ycolor','w') 
%         set(gca,'YTickLabel',[])
%         set(gca,'YTick',[])
        hold off
        nexttile(6)
        set(gca,'Color', [0 0 0])
        plot(xAxVNR, stdv, 'Color', [.5 .5 .5])
        set(gca,'Color', [0 0 0])
        hold on
        xline(xax(find(xax>3,1) + j), 'w')
        xlim([3 4])
        xlabel('time (s)')
        ylabel('VNR')
        set(gca,'box', 'off')
        set(gca,'xcolor','w')
        set(gca,'ycolor','w') 
        set(gca,'YTickLabel',[])
        set(gca,'YTick',[])
         hold off
        drawnow
       
        MM(j) = getframe(gcf);
    end
%     v = VideoWriter(['X:\Lab\Labmembers\Urs Boehm\Data\7-zArchon\composite_movie_color.avi'], 'Uncompressed AVI');
%     open(v)
%     writeVideo(v,MM)
%     close(v)
%% prepare data for fig 2 CD
% sel = V2as(ismember(V2as, find(logSR>0)));
selv3 = v3(ismember(v3, find(logSR>0)));
sel = find(ismember([data.cellId], [data(selv3).cellId]));
v3Voltage = cat(2, data(sel).trace);
v3Voltage = timeseries(v3Voltage, xAxImg(500:size(tr,1)+499));
% v3Voltage = idealfilter(v3Voltage, [2 200], 'pass');
trials = unique([data(sel).stdV]);

%% manaually defined bursts with following strong inhibition, plot traces for each cell
inhBouts = {[53 70 117 114 160 166],[],[13 122 136 183 215],[]...
    [59 72 83 94 114 124 126 128 138 139],[],[47 55 82],[]...
    [27 39 65 99 114 123 135 152 181],[],[41 94 113 139 151],[]...
    [95],[],[86 94 119 139 146 167 174 188 231],[]...
    [],[],[],[]...
    [66 78 121 129 147 152 157 160 167],[],[73 93 119 141 151 174 183 204 227 ],[]...
    [107 126 136],[],[90 102 116 129 159],[]...
    [],[],[],[]...
    [72 81 135],[],[31 64 ],[]...
    [113 129 164 185 195],[],[53 76 95 159 201 212],[]...
    [82 98 108 152 155],[],[53 104 194],[]};

selectedCells = [7 8 11 15 21 24 26 31 39 41 48 50 77 79 80 81 85];
selectedCellsAv = [];
selectedCellsVNR = [];
selectedCellsEndVNR = [];
selectedCellsEndAv = [];
curCell = [];
fN = 99;
allData = struct();
iteriter = 1;
allVNR = {};
bTs = [];
% iterate over FOVs
for everything = 1:4:11*4
    everyY = {};
    everyEndY = {};
    everyYv = [];
    everyEndYv = [];
    ind = 1;
    iter = 1;
    inhiter = [everything everything + 2];
    % iterate over trials 1 & 3 in every FOV
    for k = trials(inhiter)
        inh = inhBouts{inhiter(iter)};
        if ~isempty(inh)
            % get index for data
            tmpsel = sel(ismember([data(sel).stdV], k));
            % get index for v3Voltage
            tmpind = find(ismember([data(sel).stdV], k));
    
            % set up removed spike variable
            trRemS = zeros(size(v3Voltage.Data,1), length(tmpsel));
            trSp = zeros(size(trRemS));
            rawTrace = [];
            % iterate over cells
            for n = 1:length(tmpsel)
                % get spikes
                tmpS = data(tmpsel(n)).spikes;
                % get voltage trace
                tmpTr = (real(v3Voltage.Data(:,tmpind(n))));
                % initialize removed spike indicator
                tmpSpR = zeros(size(tmpTr));
                % keep unmodified voltage trace
                rawTrace = cat(2, rawTrace, tmpTr);
                % remove spikes
                for s = 1:length(tmpS)
                    p1 = max([tmpS(s)-3 1]);
                    p2 = min([tmpS(s)+3 length(tmpTr )]);
                    for p = p1:p2
                        tmpTr(p) = median(tmpTr(max([1 p-5]):p));
                    end
                    tmpSpR(p1:p2) = 1;
                end
                trRemS(:,n) = tmpTr;
                trSp(:,n) = tmpSpR;
            end
            %get burst in imaging time
            bE = data(tmpsel(1)).burstImg(data(tmpsel(1)).boutEndIdx)-500;
            % get burst end in VNR time
            bEVNR = data(tmpsel(1)).burstVNR(data(tmpsel(1)).boutEndIdx);
            % get all bursts and correct timing
            bT = data(tmpsel(1)).burstImg-500;
            % remove bursts that happened before imaging starts
            bT(bT<1) = [];
            bT(bT>length(v3Voltage.Time)) = [];
            
            allMean = [];
            tmpVNR =stdVAll{k};
            srate = length(tmpVNR)/25;
            xAxVNR = 1/srate:1/srate:length(tmpVNR)/srate;
            tmin = 300;
            tmax = 300;
            % iterate over cells
            for n = 1:length(tmpsel)
                currCell = data(tmpsel(n)).cellId;
                lrPos = data(tmpsel(n)).yzproj(:,1);
                allyy = [];
                allEndy = [];
                allyv = [];
                allEndyv = [];
                tmpTr = trRemS(:,n);
                tmpSpR = trSp(:,n);
                currEnd = 0;

                % iterate over inhibition events
                for p = 1:length(inh)
                    if n == 1
                        tmpbts = nan(11,1);
                        goleft = max([inh(p)-5 1]) - inh(p);
                        goright = min([inh(p)+5 length(bT)]) -inh(p);
                        tmpbts(6+goleft:6+goright) = bT(inh(p)+goleft:inh(p)+goright);
                        bTs = [bTs tmpbts];
                    end
                    yy = tmpTr(bT(inh(p))-tmin:bT(inh(p))+tmax);
                    yys = tmpSpR(bT(inh(p))-tmin:bT(inh(p))+tmax);
                    xx = v3Voltage.Time(bT(inh(p))-tmin:bT(inh(p))+tmax);
                    xx = xx-xx(tmin);
                    yy = yy-mean(yy);
                    nextEnd = bE(find((bE - bT(inh(p)))>0,1));

                    if nextEnd ~= currEnd
                        endy = tmpTr(nextEnd-tmin:nextEnd+tmax);
                        endy = endy-mean(endy);
                        allEndy = [allEndy endy];
                        currEnd = nextEnd;
                        nextEndVNR = bEVNR(find((bEVNR - data(tmpsel(1)).burstVNR(inh(p)))>0,1));
                        endyv = tmpVNR(nextEndVNR-tmin*srate/1000:nextEndVNR+tmax*srate/1000);
                        allEndyv = [allEndyv endyv];
                    end

                    allyy = [allyy yy];
                    
                    
                    yv = tmpVNR(data(tmpsel(1)).burstVNR(inh(p))-tmin*srate/1000:data(tmpsel(1)).burstVNR(inh(p))+tmax*srate/1000);
                    xv = xAxVNR(data(tmpsel(1)).burstVNR(inh(p))-tmin*srate/1000:data(tmpsel(1)).burstVNR(inh(p))+tmax*srate/1000);
                    xv = xv-xv(tmin*srate/1000);

                    allyv = [allyv yv];
                    allData(iteriter).yy = yy;
                    allData(iteriter).yys = yys;
                    allData(iteriter).yyRaw = real(rawTrace((bT(inh(p))-tmin:bT(inh(p))+tmax),n));
                    allData(iteriter).yv = yv;
                    allData(iteriter).xx = xx;
                    allData(iteriter).xv = xv';
                    allData(iteriter).currCell = currCell;
                    allData(iteriter).inhBout = inh(p);
                    allData(iteriter).trial = k;
                    allData(iteriter).lrPos = lrPos;
                    iteriter = iteriter + 1;
                end
                everyY{ind} = allyy;
                everyEndY{ind} = allEndy;
                ind = ind + 1;

                allMean = [allMean mean(allyy,2)];

                fN = fN + 1;
            end
           
            everyYv = [everyYv allyv];
            everyEndYv = [everyEndYv allEndyv];

            iter = iter + 1;
        end
    end
    if ~isempty(everyY)

        for n = 1:length(tmpsel)

            tmpMean = nanmean([everyY{n} everyY{length(tmpsel)+n}], 2);
            tmpEndMean = nanmean([everyEndY{n} everyEndY{length(tmpsel)+n}], 2);

              selectedCellsAv = cat(2, selectedCellsAv, tmpMean);
              selectedCellsEndAv = cat(2, selectedCellsEndAv, tmpEndMean);
              if length(mean(everyYv, 2)) == (tmin + tmax)*10 +1
              selectedCellsVNR = cat(2, selectedCellsVNR, interp1(1:(tmin + tmax)*10 +1, mean(everyYv, 2), 1:.2:(tmin + tmax)*10 +1)');
              selectedCellsEndVNR = cat(2, selectedCellsEndVNR, interp1(1:(tmin + tmax)*10 +1, mean(everyEndYv, 2), 1:.2:(tmin + tmax)*10 +1)');
              else
                  selectedCellsVNR = cat(2, selectedCellsVNR, mean(everyYv, 2));
                  selectedCellsEndVNR = cat(2, selectedCellsEndVNR, mean(everyEndYv, 2));
              end
        end
    end
    allVNR{end +1} = everyYv;
end


inhIdx = cat(1, data(selv3(selectedCells)).yzproj);
C = cmocean('balance', 100);
C = [C(71:100,:); C(1:30,:)];
inhIdx = inhIdx(:,1);
inhIdx = round(inhIdx.*100)-20;
% inhIdx(inhIdx>=.5) = 2;
% inhIdx(inhIdx<.5) = 1;
% C = [1 0 0; 0 0 1];

tmpTrace = nrm(selectedCellsAv(:,selectedCells));
figure
tiledlayout(2,1)
nexttile
hold on
% plot(xx, selectedCellsAv(:,selectedCells)-mean(selectedCellsAv(1:10,selectedCells),1), 'k')
for n = 1:length(inhIdx)
    plot(xx, tmpTrace(:,n), 'Color', C(inhIdx(n),:))
end

% plot(xv(1):1/50000:xv(end), selectedCellsVNR(:,selectedCells).*10-.4, 'k')
% stackplot(nrm(selectedCellsVNR(:,selectedCells))-1, linspace(xv(1), xv(end), length(selectedCellsVNR(:,selectedCells))))
plot(linspace(xv(1), xv(end), length(selectedCellsVNR(:,selectedCells))), nrm(selectedCellsVNR(:,selectedCells))-1, 'k')

nexttile
% plot(xx, selectedCellsEndAv(:,selectedCells)-mean(selectedCellsEndAv(1:10,selectedCells),1), 'k')
plot(xx, nrm(selectedCellsEndAv(:,selectedCells)), 'k')
hold on
% plot(xv(1):1/50000:xv(end), selectedCellsEndVNR(:,selectedCells).*10-.4, 'k')
plot(linspace(xv(1), xv(end), length(selectedCellsVNR(:,selectedCells))), nrm(selectedCellsEndVNR(:,selectedCells))-1, 'k')
figure
fr = 1./diff(bTs);
fr = fr.*1000;
fr(fr<15) = NaN;
plot(-4:5,fr)
hold on
plot(-4:5, mean(fr,2,'omitnan'), 'k', 'Linewidth', 2)
ylabel('cycle frequency (Hz)')
xlabel('cycle number')
%% fig 2CD
cellids = unique([allData.currCell]);
nCells = length(cellids);
% figure
% tiledlayout(ceil(sqrt(nCells)), ceil(sqrt(nCells)), 'Padding', 'none', 'Tilespacing', 'none')

for n = 1:nCells
    %     nexttile
    if ismember(n, selectedCells(9))
        trialid = [allData([allData.currCell] == cellids(n)).trial];
        trialid = trialid(1);
        boutid = [allData([allData.currCell] == cellids(n)).inhBout];
        boutid = boutid(9);%9
%         examplecid = [allData([allData.inhBout] == boutid & [allData.trial] == trialid & ismember([allData.currCell], cellids(selectedCells))).currCell];
        examplecid = [allData([allData.inhBout] == boutid & [allData.trial] == trialid).currCell];
        
%         fluo = [allData([allData.inhBout] == boutid & [allData.trial] == trialid & ismember([allData.currCell], cellids(selectedCells))).yyRaw]; 
        fluo = [allData([allData.inhBout] == boutid & [allData.trial] == trialid ).yyRaw]; 
%         fluotime = [allData([allData.inhBout] == boutid & [allData.trial] == trialid & ismember([allData.currCell], cellids(selectedCells))).xx];
        fluotime = [allData([allData.inhBout] == boutid & [allData.trial] == trialid).xx];
%         swimming = [allData([allData.inhBout] == boutid & [allData.trial] == trialid & ismember([allData.currCell], cellids(selectedCells))).yv];
        swimming = [allData([allData.inhBout] == boutid & [allData.trial] == trialid).yv];
%         swimmingtime = [allData([allData.inhBout] == boutid & [allData.trial] == trialid & ismember([allData.currCell], cellids(selectedCells))).xv];
        swimmingtime = [allData([allData.inhBout] == boutid & [allData.trial] == trialid).xv];
        figure
        tiledlayout(3,1)
        nexttile
        hold on
        c1 = 20;%16 20
        c2 = 14;%14
        yzp = data([data.cellId] == examplecid(c1)).yzproj;
        yzp = round(yzp(1)*100)-20;
        plot(fluotime(:,c1), fluo(:,c1), 'Color', C(yzp,:))
        yzp = data([data.cellId] == examplecid(c2)).yzproj;
        yzp = round(yzp(1)*100)-20;
        plot(fluotime(:,c2), fluo(:,c2), 'Color', C(yzp,:))
        plot(swimmingtime(:,1), swimming(:,1).*100-3, 'k')
        xline(0)
        xlim([-.3 .3])
        nexttile
        %         hold on
        %         fluo = [allData([allData.currCell] == cellids(n)).yy];
        %         plot(fluotime, fluo, 'Color', [.5 .5 .5])
        %         plot(fluotime, mean(fluo, 2), 'k', 'LineWidth', 2)
        %         plot(swimmingtime, swimming.*100-3, 'k')
        %         xlim([-.02 .06])
        lrIdx = cat(1, data(selv3(selectedCells)).yzproj);
        lrIdx = lrIdx(:,1);
        tmpTrace2 = nrm(selectedCellsAv(:,selectedCells));
        tmpTrace2 = tmpTrace2(300:400,:);
        inhIdx2 = zeros(size(tmpTrace2,2), 1);
        for k = 1:size(tmpTrace2,2)
            inhIdx2(k) = find(tmpTrace2(:,k)<.5, 1);
        end
       
        histogram(xx(inhIdx2)-xx(1), [0:0.003:0.035]);
        a1 = gca;
        nexttile
        hold on
        for k = 1:length(inhIdx)
            plot(xx, tmpTrace(:,k), 'Color', C(inhIdx(k),:))
        end
        
        plot(linspace(xv(1), xv(end), length(selectedCellsVNR(:,selectedCells))), nrm(selectedCellsVNR(:,selectedCells))-1, 'k')
        a2 = gca;
        linkaxes([a1 a2], 'x')
        xlim([-.02 .06])
        ylim([-1.1 1.1])
        %         hold on
        %         stackplot(fluo, fluotime)
        %         stackplot(swimming, swimmingtime)
        %         xline(0)
    end
end
%% prepare for fig s5
selector = false(size(data));
for n = 1:length(data)
    selector(n) = data(n).trial==1 &&...
        data(n).fishId==4 &&...
        data(n).zPos==1 &&...
        data(n).xypos==1 &&...
        ~isempty(data(n).av);
end

cId1 = cat(1, data(selector).cellId);

selector = false(size(data));
for n = 1:length(data)
    selector(n) = data(n).trial==1 &&...
        data(n).fishId==4 &&...
        data(n).zPos==1 &&...
        data(n).xypos==1 &&...
        ~isempty(data(n).av);
end

cId2 = cat(1, data(selector).cellId);

cId = ismember(cId1, cId2);

cId = cId1(cId);

selector = false(size(data));
for n = 1:length(data)
    selector(n) = (data(n).trial==1 ||...
        data(n).trial==1) &&...
        ismember(data(n).cellId, cId);
end


av = cat(2,data(selector).av);
amp = max(av) - min(av);
ids = find(selector);
% ids = ids(amp>=.5);

tr = cat(2,data(ids).trace);
xPos = cat(1,data(ids).xyzpos);
xPos = xPos(:,1);
yPos = cat(1,data(ids).yzproj);
yPos = yPos(:,1);
spikes = cat(1,{data(ids).spikes});
%%
% figure
% hold on
trRemS = tr;
% trRemS = tmpTr;
for n = 1:size(tr,2)
        tmpS = spikes{n};
    for k = 1:length(spikes{n})
        try
        for p = -3:3
            trRemS(tmpS(k)+p,n) = median(trRemS(max([1 tmpS(k)+p-5]):min([tmpS(k)+p size(tr,1)]),n));
        end
        catch
        end
    end
%     plot(tr(:,n)+n*5,'k')
%     plot(trRemS(:,n)+5*n,'r')
end
% figure
% stackplot(trRemS)
% hold on
% stackplot(tr)
%%
[~, posIdx] = sort(yPos);
traces = trRemS;
% traces = medfilt1(trRemS,5);
% traces = medfilt1(trRemS(:,posIdx),5);
% traces = traces - medfilt1(traces, 2001);
% traces = traces(100:end-5000,:);
% traces = medfilt1(tr,5);
% traces = tr;
[coeff, score, latent, ~, explained] = pca(traces);

%% fig. s5
srImg = 996.3;
xAxImg = 1/srImg:1/srImg:25000/srImg;
[coeff, score, latent, ~, explained] = pca(zscore(trRemS));
figure
tiledlayout(2,1)
nexttile
stackplot(score(1:10:end,1:5), xAxImg(501:10:501+length(tr)-1))
f1 = gcf;
f1.Renderer = 'painters';
nexttile
plot(xAxImg(501:10:501+length(tr)-1),trRemS(1:10:end,:))
hold on
plot(xAxImg(501:10:501+length(tr)-1), mean(trRemS(1:10:end,:),2), 'k', 'LineWidth', 2)
figure
plot(explained, 'k')
ylim([0 100])
xlim([1 22])