%% load data
opts = spreadsheetImportOptions("NumVariables", 4);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:D71";

% Specify column names and types
opts.VariableNames = ["file", "protocol", "genotype", "NMDA"];
opts.VariableTypes = ["string", "categorical", "categorical", "double"];

% Specify variable properties
opts = setvaropts(opts, "file", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["file", "protocol", "genotype"], "EmptyFieldRule", "auto");

% Import the data
master = readtable("..\Data\ChrMaster.xlsx", opts, "UseExcel", false);
clear opts
load('..\Data\NMDAincreasingConc.mat');
%% split complete protocols
datas2 = struct();
master2 = table();
idx = 1;
for n = 1:length(datas)
    if master.protocol(n) == 'complete'
        T = datas(n).table;
        time = T.time - T.time(1);
        T1 = T(time<=360,:);
        T2 = T(time>360,:);
        datas2(idx).table = T1;
        master2(idx,:) = master(n,:);
        master2.protocol(idx) = 'const';
        master2.NMDA(idx) = master.NMDA(n);
        idx = idx + 1;
        datas2(idx).table = T2;
        master2(idx,:) = master(n,:);
        master2.protocol(idx) = 'triggered';
        master2.NMDA(idx) = master.NMDA(n);
        idx = idx + 1;
    else
        datas2(idx).table = datas(n).table;
        master2(idx,:) = master(n,:);
        idx = idx + 1;
    end
end
%% assign fish id
for n = 1:height(master2)
    id = str2num(regexp(master2.file(n), '202[0-1][0-9]*', 'match'))*10;
    id = id + str2num(regexp(master2.file(n), '(?<=\\l)[0-9]*','match'));
    master2.id(n) = id;
end
%% analyze swimming
id = 1;
TT = struct();
for i =  1:length(datas2)
    data = datas2(i).table;
    %make structure for later
    T=struct('Bout_starts',[],'Bout_ends',[],'Bout_duration',[],'Interbout_starts',[],...
        'Interbout_ends',[],'Interbout_duration',[],'N_cycles',[],'TBF',[],'Max_amp',[],...
        'Mean_amp',[], 'Laterality',[], 'allPeaks', []);
    ang = data.angle-mean(data.angle);
    ang = lowpass(ang, 25, 200);
    ang = ang-medfilt1(ang, 250);
    
    %finds local minima/maxima
    [lmax,lmin]=peakdet(ang, 0.03); % ADJUST second value to change sensitivity
    
    
    %Set parameters for bouts and interbouts, this part and below modified from fictive code bout detection!
    interbout_Int=24;
    %interbout_Int=15;  %minimum number of frames between bouts  (15=60ms)
    min_bout_cutoff=7; %minimum length of bout  (7=28ms)
    
    %Here combine lmax and lmin to make a list of all detected peaks so doesn't exclude one side of zero
    %Also add the last detected peak as last frame so doesn't exclude last bout later.
    if ~isempty(lmax) && ~isempty(lmin)
        allpeaks=union(lmax(:,1),lmin(:,1));
    elseif isempty(lmax) && ~isempty(lmin)
        allpeaks = lmin(:,1);
    elseif ~isempty(lmax) && isempty(lmin)
        allpeaks = lmax(:,1);
    else
        allpeaks = [];
    end
    if size(allpeaks,1)==1
        allpeaks = allpeaks';
    end
    instFreq = (1./(diff(allpeaks)));
    instFreq = [0; instFreq];
    highFreqFlag = instFreq>0.5;
    allpeaks(highFreqFlag) = [];
    instFreq = (1./(diff(allpeaks)));
    instFreq = [0; instFreq];
    
    Diff_thresh=1./instFreq;
    Find_Interbouts=find(Diff_thresh>interbout_Int);
    %IMPORTANT: Here, Find_Interbouts is a list of Diff_thresh indexes which contain the differences between sequential timepoints, not a list of timepoints itself
    %Remember, the first high difference value (ie Find_Interbouts(1) is going to be the last thresholded point on the first bout!
    if ~isempty(allpeaks)
        T.Bout_starts = allpeaks(Find_Interbouts);
        T.Bout_ends = [allpeaks(Find_Interbouts(2:end)-1); allpeaks(end)];
        for bouts = 1:length(T.Bout_starts)
            bs = find(ang(1:T.Bout_starts(bouts))<std(ang) & ang(1:T.Bout_starts(bouts))>-std(ang),1, 'last');
            if ~isempty(bs)
                T.Bout_starts(bouts) = bs;
            end
            be = find(ang(T.Bout_ends(bouts):end)<std(ang) & ang(T.Bout_ends(bouts):end)>-std(ang),1, 'first');
            if ~isempty(be)
                be = be + T.Bout_ends(bouts) - 1;
                T.Bout_ends(bouts) = be;
            end
        end
        T.Bout_duration = T.Bout_ends-T.Bout_starts;
        
        T.Bout_starts(T.Bout_duration<=min_bout_cutoff) = [];
        T.Bout_ends(T.Bout_duration<=min_bout_cutoff) = [];
        
        T.Bout_duration(T.Bout_duration<=min_bout_cutoff) = [];
        
        T.allPeaks = allpeaks;
    end
    framerate=200;
    
    T.Bout_duration=(T.Bout_duration./framerate)*1000;
    
    %Calculate things about each bout
    %Ncycles, Avg TBF, Max Amp? Laterality?
    %Laterality:  0=Right, 1=Left
    
    for j=1:length(T.Bout_starts)
        boutpeaks=find(allpeaks>=T.Bout_starts(j)&allpeaks<=T.Bout_ends(j));
        T.N_cycles(j)=length(boutpeaks)/2;
        T.TBF(j)=1000*T.N_cycles(j)/T.Bout_duration(j);
        T.Max_amp(j)=max(abs(ang(allpeaks(boutpeaks))));
        T.vigor(j) = mean(data.vigor(T.Bout_starts(j):T.Bout_ends(j)));
        T.vigor_sum(j) = sum(data.vigor(T.Bout_starts(j):T.Bout_ends(j)));
        if j<length(T.Bout_starts)
            T.vigor_tot(j) = mean(data.vigor(T.Bout_starts(j):T.Bout_starts(j+1)));
        else
            T.vigor_tot(j) = mean(data.vigor(T.Bout_starts(j):end));
        end
        if length(boutpeaks) >1
            T.Mean_amp(j) = mean([abs(ang(allpeaks(boutpeaks(1)))); abs(diff(ang(allpeaks(boutpeaks))))./2]);
        else
            T.Mean_amp(j) = abs(ang(allpeaks(boutpeaks)));
        end
        if ang(allpeaks(boutpeaks(1)))>0
            T.Laterality(j)=-1;
        else
            T.Laterality(j)=1;
        end
        clear boutpeaks
    end
    TT(id).T = T;
    id = id + 1;
end

sel = 1:length(TT);
for i = 1:length(TT)
    T = TT(i).T;
    data = datas2(sel(i)).table;
    ang = data.angle-mean(data.angle);
    allpeaks = T.allPeaks;
    ags = {};
    pks = {};
    ls = {};
    lsm = [];
    OMR = [];
    gain = [];
    stdp = [];
    intensity = [];
    xAx = 1/200:1/200:length(ang)/200;
    laser = data.laser;
    for n = 1:length(T.Bout_starts)
        tmp2 = T.Laterality(n).*ang(T.Bout_starts(n):min([T.Bout_ends(n) length(ang)]));
        ags{n} = tmp2;
        pks{n} = allpeaks(allpeaks>=T.Bout_starts(n) & allpeaks<=(T.Bout_ends(n))) - T.Bout_starts(n) +1;
        stdp(n) = std(sqrt(ags{n}(pks{n}).^2));
        ls{n} = laser(T.Bout_starts(n):T.Bout_ends(n));
        if max(ls{n}) > 0
            ls{n} = ls{n}./max(ls{n});
        end
        lsm(n) = mean(ls{n});
        OMR(n) = mean(data.OMR(T.Bout_starts(n):T.Bout_ends(n)));
        gain(n) = mean(data.gain(T.Bout_starts(n):T.Bout_ends(n)));
        intensity(n) = max(data.laser(T.Bout_starts(n):T.Bout_ends(n)));
    end
    T.laser_ratio = lsm;
    lsm(lsm>=0.5) = 1;
    lsm(lsm<0.5) = 0;
    intensity = intensity.*lsm;
    T.ags = ags;
    T.pks = pks;
    T.lsm = lsm;
    T.OMR = OMR;
    T.gain = gain;
    T.stdp = stdp;
    T.intensity = intensity;
    TT(i).T = T;
end
%% generate long table
sw = 0;
for n = 1:length(TT)
    sw = sw+length(TT(n).T.Bout_starts);
end
longTable = table('Size', [sw, 16], 'variableTypes', {'double', 'double',...
    'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double',...
    'categorical', 'double', 'categorical'}, 'variableNames', {'Duration', 'Ncycles', 'TBF',...
    'vigor', 'vigorTot', 'fictiveDistance', 'laserRatio', 'MaxAmp', 'MeanAmp', 'SumVigor', 'fishID', 'NMDAconc', 'Laser', 'Protocol', 'Intensity', 'Genotype'});
id = 1;
for n = 1:length(TT)
    for k = 1:length(TT(n).T.Bout_starts)
        longTable.Duration(id) = TT(n).T.Bout_duration(k);
        longTable.Ncycles(id) = TT(n).T.N_cycles(k);
        longTable.TBF(id) = TT(n).T.TBF(k);
        longTable.MaxAmp(id) = TT(n).T.Max_amp(k);
        longTable.MeanAmp(id) = TT(n).T.Mean_amp(k);
        longTable.fishID(id) = master2.id(n);
        longTable.Genotype(id) = master2.genotype(n);
        longTable.Laser(id) = TT(n).T.lsm(k);
        longTable.Protocol(id) = master2.protocol(n);
        longTable.Intensity(id) = TT(n).T.intensity(k);
        longTable.NMDAconc(id) = master2.NMDA(n);
        longTable.vigor(id) = TT(n).T.vigor(k);
        longTable.vigorTot(id) = TT(n).T.vigor_tot(k);
        longTable.fictiveDistance(id) = TT(n).T.vigor_sum(k);
        longTable.laserRatio(id) = TT(n).T.laser_ratio(k);
        id = id + 1;
    end
end
%% larvae to exclude
exclude = [202101124, 202101123, 202012228, 202012227, 202012221];

%% scatters per fish
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% individual NMDA conc for each fish
% fish id   NMDA conc
indvConc =  [202101126 125;...
202101125 125;...
202012229 100;...
202012230 100;...% (little swimming at led off)
202101122 100;...
202101121 125;...% (little swimming at led off)
212101113 50;...%  (little swimming at led off, at higher NMDA probability goes >.8)
202012226 100;...% (little swimming at led off)
202101112 75;...
202101111 75;...
202012225 100];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
mnrmMaOff = [];
mnrmMaOn = [];
valsAll = [];
laserAll = [];
fidAll = [];
id = 1;
for n = unique(longTable.fishID)'
    conc = indvConc(indvConc(:,1)==n,2);
    if isempty(conc)
        conc = 0;
    end
    valsOff = log10(longTable.fictiveDistance(longTable.Laser==0 & longTable.Ncycles>.5 & longTable.NMDAconc==100 & longTable.fishID==n & ~ismember(longTable.fishID, exclude) & ~(longTable.laserRatio>.3 & longTable.laserRatio<.7)));
    valsOn = log10(longTable.fictiveDistance(longTable.Laser==1 & longTable.Ncycles>.5 & longTable.NMDAconc==100 & longTable.fishID==n & ~ismember(longTable.fishID, exclude) & ~(longTable.laserRatio>.3 & longTable.laserRatio<.7)));
%     valsOff = nrmMa(sel & longTable.Laser==0 & longTable.fishID==n);
%     valsOn = nrmMa(sel & longTable.Laser==1 & longTable.fishID==n);
    if length(valsOff)>5
%         disp(num2str(n))
        mnrmMaOff(id) = mean(valsOff);
        mnrmMaOn(id) = mean(valsOn);
        valsAll = [valsAll; valsOn; valsOff];
        laserAll = [laserAll; ones(size(valsOn)); zeros(size(valsOff))];
        fidAll = [fidAll; ones(size(valsOn)).*n; ones(size(valsOff)).*n];
%         disp(['On: ' num2str(length(valsOn)) ' Off: ' num2str(length(valsOff))])
        id = id + 1;
    end
end

fittable = table(valsAll, laserAll, fidAll, 'VariableNames',...
    {'fictiveDistance', 'Laser', 'fishID'});
lmeDist = fitlme(fittable, 'fictiveDistance ~ Laser + (Laser|fishID)');

figure(3)
subplot(1,4,1)
ph = plotSpread(([(mnrmMaOff) (mnrmMaOn)]), 'distributionIdx', [ones(size(mnrmMaOff)) ones(size(mnrmMaOn)).*2],...
    'distributionColor', {'k', 'b'}, 'distributionMarkers', 'o', 'xNames', {'LED Off', 'LED On'});
% title('mean fictive Distance at 100 \muM NMDA')
ylabel('log_{10}(mean fictivce Distance) (a.u.)')
% ylim([0 max(ylim(ph{3}))])
ph = get(ph{3}, 'Children');
for n = 1:length(ph)
    if strcmp(ph(n).DisplayName, '1')
        ph(n).MarkerFaceColor = 'k';
    elseif strcmp(ph(n).DisplayName, '2')
        ph(n).MarkerFaceColor = 'b';
    end
end
plot(([mnrmMaOff; mnrmMaOn]), 'k')

% figure(4)
% s = subplot(1,4,1);
% scatter(mnrmMaOff, mnrmMaOn, 'k', 'filled')
% hold on
% line([min([s.YLim s.XLim]) max([s.YLim s.XLim])], [min([s.YLim s.XLim]) max([s.YLim s.XLim])], 'Color', 'k')
% daspect([1 1 1])

mnrmDurOff = [];
mnrmDurOn = [];
selId = [];
valsAll = [];
laserAll = [];
fidAll = [];
id = 1;
for n = unique(longTable.fishID)'
    conc = indvConc(indvConc(:,1)==n,2);
    if isempty(conc)
        conc = 0;
    end
    valsOff = (longTable.Duration(longTable.Laser==0 & longTable.Ncycles>.5 & longTable.NMDAconc==100 & longTable.fishID==n & ~ismember(longTable.fishID, exclude)));
    valsOn = (longTable.Duration(longTable.Laser==1 & longTable.Ncycles>.5 & longTable.NMDAconc==100 & longTable.fishID==n & ~ismember(longTable.fishID, exclude)));
    if length(valsOff)>5
        mnrmDurOff(id) = mean(valsOff);
        mnrmDurOn(id) = mean(valsOn);
        selId(id) = n;
        valsAll = [valsAll; valsOn; valsOff];
        laserAll = [laserAll; ones(size(valsOn)); zeros(size(valsOff))];
        fidAll = [fidAll; ones(size(valsOn)).*n; ones(size(valsOff)).*n];
        id = id + 1;
    end
end

fittable = table(log10(valsAll), laserAll, fidAll, 'VariableNames',...
    {'Duration', 'Laser', 'fishID'});
lmeDur = fitlme(fittable, 'Duration ~ Laser + (Laser|fishID)');

figure(3)
subplot(1,4,2)
ph = plotSpread([(mnrmDurOff) (mnrmDurOn)], 'distributionIdx', [ones(size(mnrmDurOff)) ones(size(mnrmDurOn)).*2],...
    'distributionColor', {'k', 'b'}, 'distributionMarkers', 'o', 'xNames', {'LED Off', 'LED On'});
title('mean bout duration')
ylabel('time (ms)')
% ylim([0 max(ylim(ph{3}))])
ph = get(ph{3}, 'Children');
for n = 1:length(ph)
    if strcmp(ph(n).DisplayName, '1')
        ph(n).MarkerFaceColor = 'k';
    elseif strcmp(ph(n).DisplayName, '2')
        ph(n).MarkerFaceColor = 'b';
    end
end
plot([mnrmDurOff; mnrmDurOn], 'k')

% figure(4)
% s = subplot(1,4,2);
% scatter(mnrmDurOff, mnrmDurOn, 'k', 'filled')
% hold on
% line([min([s.YLim s.XLim]) max([s.YLim s.XLim])], [min([s.YLim s.XLim]) max([s.YLim s.XLim])], 'Color', 'k')
% daspect([1 1 1])

mnrmTBFOff = [];
mnrmTBFOn = [];
valsAll = [];
laserAll = [];
fidAll = [];
id = 1;
for n = unique(longTable.fishID)'
    conc = indvConc(indvConc(:,1)==n,2);
    if isempty(conc)
        conc = 0;
    end
    valsOff = longTable.TBF(longTable.Laser==0 & longTable.Ncycles>.5 & longTable.NMDAconc==100 & longTable.fishID==n & ~ismember(longTable.fishID, exclude) & ~(longTable.laserRatio>.3 & longTable.laserRatio<.7));
    valsOn = longTable.TBF(longTable.Laser==1 & longTable.Ncycles>.5 & longTable.NMDAconc==100 & longTable.fishID==n & ~ismember(longTable.fishID, exclude) & ~(longTable.laserRatio>.3 & longTable.laserRatio<.7));
    if length(valsOff)>5
        mnrmTBFOff(id) = mean(valsOff);
        mnrmTBFOn(id) = mean(valsOn);
        valsAll = [valsAll; valsOn; valsOff];
        laserAll = [laserAll; ones(size(valsOn)); zeros(size(valsOff))];
        fidAll = [fidAll; ones(size(valsOn)).*n; ones(size(valsOff)).*n];
        id = id + 1;
    end
end
fittable = table(valsAll, laserAll, fidAll, 'VariableNames',...
    {'TBF', 'Laser', 'fishID'});
lmeTBF = fitlme(fittable, 'TBF ~ Laser + (Laser|fishID)');

figure(3)
subplot(1,4,3)
ph = plotSpread([mnrmTBFOff mnrmTBFOn], 'distributionIdx', [ones(size(mnrmTBFOff)) ones(size(mnrmTBFOn)).*2],...
    'distributionColor', {'k', 'b'}, 'distributionMarkers', 'o', 'xNames', {'LED Off', 'LED On'});
title('mean beat frequency')
ylabel('frequency [a.u.]')
% ylim([0 max(ylim(ph{3}))])
ph = get(ph{3}, 'Children');
for n = 1:length(ph)
    if strcmp(ph(n).DisplayName, '1')
        ph(n).MarkerFaceColor = 'k';
    elseif strcmp(ph(n).DisplayName, '2')
        ph(n).MarkerFaceColor = 'b';
    end
end
plot([mnrmTBFOff; mnrmTBFOn], 'k')
ylim([0 19])
% figure(4)
% s = subplot(1,4,3);
% scatter(mnrmTBFOff, mnrmTBFOn, 'k', 'filled')
% hold on
% line([min([s.YLim s.XLim]) max([s.YLim s.XLim])], [min([s.YLim s.XLim]) max([s.YLim s.XLim])], 'Color', 'k')
% daspect([1 1 1])


mnrmAmpOff = [];
mnrmAmpOn = [];
valsAll = [];
laserAll = [];
fidAll = [];
id = 1;
for n = unique(longTable.fishID)'
    conc = indvConc(indvConc(:,1)==n,2);
    if isempty(conc)
        conc = 0;
    end
    valsOff = longTable.vigor(longTable.Laser==0 & longTable.Ncycles>.5 & longTable.NMDAconc==100 & longTable.fishID==n & ~ismember(longTable.fishID, exclude) & ~(longTable.laserRatio>.3 & longTable.laserRatio<.7));
    valsOn = longTable.vigor(longTable.Laser==1 & longTable.Ncycles>.5 & longTable.NMDAconc==100 & longTable.fishID==n & ~ismember(longTable.fishID, exclude) & ~(longTable.laserRatio>.3 & longTable.laserRatio<.7));
    if length(valsOff)>5
        mnrmAmpOff(id) = mean(valsOff);
        mnrmAmpOn(id) = mean(valsOn);
        valsAll = [valsAll; valsOn; valsOff];
        laserAll = [laserAll; ones(size(valsOn)); zeros(size(valsOff))];
        fidAll = [fidAll; ones(size(valsOn)).*n; ones(size(valsOff)).*n];
        id = id + 1;
    end
end

fittable = table(valsAll, laserAll, fidAll, 'VariableNames',...
    {'Amp', 'Laser', 'fishID'});
lmeAmp = fitlme(fittable, 'Amp ~ Laser + (Laser|fishID)');

figure(3)
subplot(1,4,4)
ph = plotSpread([mnrmAmpOff mnrmAmpOn], 'distributionIdx', [ones(size(mnrmAmpOff)) ones(size(mnrmAmpOn)).*2],...
    'distributionColor', {'k', 'b'}, 'distributionMarkers', 'o', 'xNames', {'LED Off', 'LED On'});
title('mean amplitude')
ylabel('swim strength (a.u.)')
% ylim([0 max(ylim(ph{3}))])
ph = get(ph{3}, 'Children');
for n = 1:length(ph)
    if strcmp(ph(n).DisplayName, '1')
        ph(n).MarkerFaceColor = 'k';
    elseif strcmp(ph(n).DisplayName, '2')
        ph(n).MarkerFaceColor = 'b';
    end
end
plot([mnrmAmpOff; mnrmAmpOn], 'k')


ratio_dur = mean(mnrmDurOn./mnrmDurOff)
sem_dur = std(mnrmDurOn./mnrmDurOff)/sqrt(numel(mnrmDurOff))

ratio_strength = mean(mnrmAmpOn./mnrmAmpOff)
sem_strength = std(mnrmAmpOn./mnrmAmpOff)/sqrt(numel(mnrmAmpOff))

ratio_freq = mean(mnrmTBFOn./mnrmTBFOff)
sem_freq = std(mnrmTBFOn./mnrmTBFOff)/sqrt(numel(mnrmTBFOff))
%% plot repetitions on top of each other

s11 = [];
id = 1;
nmdaConc = [];
for n = 1:length(datas2)
    if ~ismember(master2.id(n), exclude)
        ang = zeros(size(datas2(n).table.vigor));
        for k = 1:length(TT(n).T.Bout_starts)
            ang(TT(n).T.Bout_starts(k):TT(n).T.Bout_ends(k)) = 1;
%             ang(TT(n).T.Bout_starts(k):TT(n).T.Bout_starts(k)+10) = 1;
%             ang(TT(n).T.Bout_ends(k):TT(n).T.Bout_ends(k)+10) = 1;
%             ang(TT(n).T.Bout_starts(k):TT(n).T.Bout_ends(k)) = datas2(n).table.vigor(TT(n).T.Bout_starts(k):TT(n).T.Bout_ends(k));
        end
%         ang = datas2(n).table.angle;
%         ang = datas2(n).table.vigor;
        laser = datas2(n).table.laser;
%         ang = ang-medfilt1(ang, 250);
%         ang = medfilt1(ang,5);
        
%         figure(master2.id(n))
%         hold on
        loff = [1; find(diff(laser)<0)];
        lon = [find(diff(laser)>0)];
        l = round(mean(diff(loff)./2));
        s1 = [];
        s2 = [];
        conc = unique(master2.NMDA);
%         subplot(4,1,find(conc == master2.NMDA(n)))
%         hold on
        for k = 1:length(lon)-1
%             plot(1/200:1/200:l*2/200, ang(lon(k)-round(l/2):lon(k)-round(l/2)+l*2-1), 'k')
%             ylim([0 0.3])
%             s1(:,k) = movmean(abs(ang(lon(k)-round(l/2):lon(k)-round(l/2)+l*2-1)), 20);
            s1(:,k) = ang(lon(k)-round(l/2):lon(k)-round(l/2)+l*2-1);
%             ylim([0 0.15])
        end
%         subplot(4,1,4)
%         xlabel('time (s)')
%         ylabel('lowpass filtered tail angle')
        
%         figure(99)
%         hold on
%         subplot(4,1,find(conc == master2.NMDA(n)))
%         hold on
%         plot(1/200:1/200:l*2/200, nanmean(s1,2), 'r')
        s11(:,id) = nanmean(s1,2);
        nmdaConc(id) = master2.NMDA(n);
        id = id + 1;
    end
end
pc = [151 135 158;  207 145 176; 245 187 171]./255;
figure
hold on
patch([1.5 4.5 4.5 1.5], [0 0 0.08 0.08], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off')
cm = cmocean('-matter', 4);
for n = 1:3
    av = nanmean(s11(:,nmdaConc==conc(n)),2);
    sem = nanstd(s11(:,nmdaConc==conc(n)),[],2)./sqrt(sum(nmdaConc==conc(n)));
%     patch([1/200:1/200:l*2/200 fliplr(1/200:1/200:l*2/200)], [av' + sem' flipud(av-sem)'],...
%         cm(n,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'Off')
    patch([1/200:1/200:l*2/200 fliplr(1/200:1/200:l*2/200)], [av' + sem' flipud(av-sem)'],...
        pc(n,:), 'EdgeColor', 'none', 'HandleVisibility', 'Off')
    plot(1/200:1/200:l*2/200, nanmean(s11(:,nmdaConc==conc(n)),2),'Color', cm(n,:), 'LineWidth', 2)
end

xlabel('time (s)')
legend({'50 \muM NMDA', '75 \muM NMDA', '100 \muM NMDA'})

%% example fig. 5C

 conc = unique(master2.NMDA);
for n = 1:length(datas2)
    if master2.id(n) == 202101112 && master2.NMDA(n) == 100
%     if master2.NMDA(n) == 100
        figure(master2.id(n))
        
        hold on
        t = datas2(n).table.time-datas2(n).table.time(n);
        ang = datas2(n).table.angle-mean(datas2(n).table.angle);
        vig = datas2(n).table.vigor;
        ang = lowpass(ang, 25, 200);
        ang = ang-medfilt1(ang, 250);
        %         ang = datas2(n).table.angle-mean(datas2(n).table.angle);
%         yyaxis left
        hold on
        patch([t; 0], [datas2(n).table.laser./10; 0]-.2, 'b', 'FaceAlpha', .1, 'EdgeColor', 'none')
        plot(t, ang, 'k')
        plot(t, vig, 'r')
        ylim([-.25 .25])
%         plot(t(TT(n).T.Bout_starts), ang(TT(n).T.Bout_starts), '.g')
%         plot(t(TT(n).T.Bout_ends), ang(TT(n).T.Bout_ends), '.r')
%         yyaxis right
%         hold on
%         bs = TT(n).T.Bout_starts;
%             be = TT(n).T.Bout_ends;
%             pks = TT(n).T.allPeaks;
%             for k = 1:length(TT(n).T.Bout_starts)
%                 
%                 pt = pks(pks>bs(k) & pks<be(k));
%                 pt = pt(1:2:end);
%                 freq = 1./(diff(pt).*0.005);
%                 plot(t(pt(2:end)), freq, '-r')
%             end
%             ylim([0 19])
        xlim([84.5 90])
    end
end