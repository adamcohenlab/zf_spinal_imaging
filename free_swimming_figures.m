% load('..\Data\OMRrepeat.mat')
%load larva info
opts = spreadsheetImportOptions("NumVariables", 2);
% Specify sheet and range
opts.Sheet = "Tabelle1";
opts.DataRange = "A2:B78";
% Specify column names and types
opts.VariableNames = ["folders", "knockout"];
opts.VariableTypes = ["string", "categorical"];
% Specify variable properties
opts = setvaropts(opts, "folders", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["folders", "knockout"], "EmptyFieldRule", "auto");
% Import the data
master = readtable("..\Data\freeSwimmingMaster.xlsx", opts, "UseExcel", false);


pxscale = 69/1325; % mm / pixel
%% calculate tail angle and kinematics from raw data
TT = struct();
for k = 1:size(rawData,2)
    behaviorlog = rawData(k).behavior;
    stimuluslog = rawData(k).stimulus;
    displacement = sqrt(behaviorlog.f0_vx.^2 + behaviorlog.f0_vy.^2);
    displint = interp1(find(~isnan(displacement)), displacement(~isnan(displacement)), 1:length(displacement));
    displint = medfilt1(displint,7);
    theta = behaviorlog.f0_theta;
    theta = interp1(find(~isnan(theta)), theta(~isnan(theta)), 1:length(theta));
    theta = medfilt1(theta,3);
    % stimvel = smooth(medfilt1([0; diff(stimuluslog.moving_gratings_x)],3));
    
    
    xpos = interp1(find(~isnan(behaviorlog.f0_x)), behaviorlog.f0_x(~isnan(behaviorlog.f0_x)), 1:length(behaviorlog.f0_x));
    ypos = interp1(find(~isnan(behaviorlog.f0_y)), behaviorlog.f0_y(~isnan(behaviorlog.f0_y)), 1:length(behaviorlog.f0_y));
    
    T=struct('Bout_starts',[],'Bout_ends',[],'Bout_duration',[],'Interbout_starts',[],...
        'Interbout_ends',[],'Interbout_duration',[],'N_cycles',[],'TBF',[],'Max_amp',[],...
        'Laterality',[], 'allPeaks', [], 'MeanSpeed', [], 'BoutStartTime', [],...
        'BoutEndTime', [], 'BoutDirection', [], 'ang', [], 'displacement', [], 'theta', [], 'omrDirection', [], 'yaw', []);
    T.displacement = displacement;
    
    T.theta = theta;
    tb = behaviorlog.t;
    ts = stimuluslog.t;
    vs = stimuluslog.conditional_velx;
    vs(isnan(vs)) = 0;
    [~,b,~] = unique(ts);
    omrInt = interp1(ts(b), vs(b), tb, 'linear', 0);
%     a = table2array(rawData(k).behavior);
%     tailAngles = a(:,8:16);
%     tailAngles(tailAngles>pi) = tailAngles(tailAngles>pi)-2*pi;
%     tailAngles(tailAngles<-pi) = tailAngles(tailAngles<-pi)+2*pi;
     ang = behaviorlog.f0_theta_07;
     angr = behaviorlog.f0_theta_00;
     ang(ang>pi) = ang(ang>pi)-2*pi;
    ang(ang<-pi) = ang(ang<-pi)+2*pi;
%     tmp = [];
%    
%     
%     for n = 2:size(tailAngles,2)
%         tmp(:,n) = 2*(tailAngles(:,n)-tailAngles(:,n-1));
%     end
%     tailAngles(:,2:end) = tmp(:,2:end);
%     
%     tailAngles = cumsum(tailAngles, 2);
%     ang = tailAngles(:,end);
       
%         ang(ang>2 | ang<-2) = NaN;
    nanflag = movmean(ang,[0 10], 'omitnan', 'Endpoints', 1);
    nankeep = false(size(nanflag));
    for nans = find(isnan(nanflag))'
        nankeep(nans:nans+10) = true;
    end
    ang(nankeep) = 0;
    ang2 = interp1(find(~isnan(ang)), ang(~isnan(ang)), 1:length(ang), 'linear', 0);
    ang = ang2-mean(ang2);
    ang = lowpass(ang, 60, 300);
%     ang = medfilt1(ang2, 3);
    vigor = movstd(ang,15);
    T.ang = ang;
    T.vigor = vigor;
%     T.ang_orig = behaviorlog.f0_theta_07;
    %finds local minima/maxima
%     [lmax,lmin]=peakdet(ang-medfilt1(ang,31), 0.25); % ADJUST second value to change sensitivity
anghp = ang-medfilt1(ang,38);
% [lmax,lmin]=peakdet(anghp, 0.25);
[~, lmax] = findpeaks(anghp, 'MinPeakHeight', 0.12, 'MinPeakProminence', 0.3, 'MinPeakDistance', 5);
[~, lmin] = findpeaks(-1.*anghp, 'MinPeakHeight', 0.12, 'MinPeakProminence', 0.3, 'MinPeakDistance', 5);

    
    %Set parameters for bouts and interbouts, this part and below modified from fictive code bout detection!
    interbout_Int=30;
    %interbout_Int=15;  %minimum number of frames between bouts  (18=60ms)
    min_bout_cutoff=6; %minimum length of bout  (6=20ms)
    
    %Here combine lmax and lmin to make a list of all detected peaks so doesn't exclude one side of zero
    %Also add the last detected peak as last frame so doesn't exclude last bout later.
    
  %Here combine lmax and lmin to make a list of all detected peaks so doesn't exclude one side of zero
    %Also add the last detected peak as last frame so doesn't exclude last bout later.
    if ~isempty(lmax) && ~isempty(lmin)
        allpeaks=union(lmax,lmin);
    elseif isempty(lmax) && ~isempty(lmin)
        allpeaks = lmin;
    elseif ~isempty(lmax) && isempty(lmin)
        allpeaks = lmax;
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
    %     
    %     lastpeak=max(allpeaks);
    %     lastindex=find(allpeaks==lastpeak);
    %     allpeaks(lastindex+1)=length(ang);
    %     

    
    Diff_thresh=1./instFreq;
    Find_Interbouts=find(Diff_thresh>interbout_Int);
    %IMPORTANT: Here, Find_Interbouts is a list of Diff_thresh indexes which contain the differences between sequential timepoints, not a list of timepoints itself
    %Remember, the first high difference value (ie Find_Interbouts(1) is going to be the last thresholded point on the first bout!
    if ~isempty(allpeaks)
        bs = allpeaks(Find_Interbouts);
        be = [allpeaks(Find_Interbouts(2:end)-1); allpeaks(end)];
        T.Bout_starts = bs(bs~=be);
        T.Bout_ends = be(bs~=be);
        
        for bouts = 1:length(T.Bout_starts)
            bs = find(sign(anghp(T.Bout_starts(bouts))).*anghp(1:T.Bout_starts(bouts))<0,1, 'last');
            if ~isempty(bs)
                T.Bout_starts(bouts) = bs;
            end
            be = find(sign(anghp(T.Bout_ends(bouts))).*anghp(T.Bout_ends(bouts):end)<0,1, 'first');
            if ~isempty(be)
                be = be + T.Bout_ends(bouts) - 1;
                T.Bout_ends(bouts) = be;
            end
        end
    end
    T.Bout_duration = T.Bout_ends-T.Bout_starts;
    T.Bout_starts(T.Bout_duration<=min_bout_cutoff) = [];
    T.Bout_ends(T.Bout_duration<=min_bout_cutoff) = [];
    
    T.Bout_duration(T.Bout_duration<=min_bout_cutoff) = [];
    
    T.allPeaks = allpeaks;
    
% figure
% plot(ang)
% hold on
% scatter(T.Bout_starts, ang(T.Bout_starts), 'g', 'filled')
% scatter(T.Bout_ends, ang(T.Bout_ends), 'r', 'filled')
% scatter(allpeaks, ang(allpeaks), 'k', 'filled')

    %
    %     %Plots bout start and ends
    %     %starts are green, ends are red
    %     ha(2)=subplot(2,1,2);
    %     plot(ang,'k');
    %     hold on;
    %     if ~isempty(T.Bout_starts)
    %     plot(T.Bout_starts, max(ang)/2, 'g.', 'Markersize', 16);
    %     hold on;
    %     plot(T.Bout_ends, max(ang)/2, 'r.', 'Markersize', 16);
    %     end
    %     plot(data.laser./10, 'b')
    %     end
    %     title('Bout Starts and Ends');
    %     linkaxes(ha,'xy');
    %     xlabel('Frames');
    %
    %Convert Durations from Frames to MS!
    framerate=300;
    
    
    T.Bout_duration=(T.Bout_duration./framerate)*1000;
    T.BoutStartTime = behaviorlog.t(T.Bout_starts);
    T.BoutEndTime = behaviorlog.t(T.Bout_ends);
    T.Interbout_duration = 1000.*[T.Bout_starts(2:end)-T.Bout_ends(1:end-1); NaN]./framerate;
    %Calculate things about each bout
    %Ncycles, Avg TBF, Max Amp? Laterality?
    %Laterality:  0=Right, 1=Left
    tmpYaw = theta-medfilt1(theta, 11);
    for j=1:length(T.Bout_starts)
        T.StimVel(j) = abs(omrInt(T.Bout_starts(j)));
        T.omrDirection(j) = sign(omrInt(T.Bout_starts(j)));
        T.Xpos_start(j) = behaviorlog.f0_x(T.Bout_starts(j));
        boutpeaks=find(allpeaks>=T.Bout_starts(j)&allpeaks<=T.Bout_ends(j));
        T.Boutpeaks{j} = allpeaks(boutpeaks);
        T.N_cycles(j)=length(boutpeaks)/2;
        T.TBF(j)=1000/(T.Bout_duration(j)/(T.N_cycles(j)));
        T.periodStd(j) = std(diff(T.Boutpeaks{j})./framerate);
        T.dperiodStd(j) = std(diff(diff(T.Boutpeaks{j})./framerate));
        T.meandPeriod(j) = mean(diff(diff(T.Boutpeaks{j})./framerate));
        T.meanPeriod(j) = mean(diff(T.Boutpeaks{j})./framerate);
%         T.Max_amp(j)=max(abs(ang(allpeaks(boutpeaks))));
        T.Max_amp(j) = max(abs(diff(ang(allpeaks(boutpeaks)))));
        T.Mean_amp(j) = mean(abs(diff(ang(allpeaks(boutpeaks)))));
        tmpBout = stimuluslog.t >= T.BoutStartTime(j) & stimuluslog.t <= T.BoutEndTime(j);
        if sum(tmpBout)==0
            tmpBout = stimuluslog.t >= T.BoutStartTime(j);
            tmpBout(find(tmpBout,1)) = 1;
        end
        T.MeanSpeed(j) = nanmean(displacement(T.Bout_starts(j):T.Bout_ends(j)));
        T.distance(j) = nansum(displacement(T.Bout_starts(j):T.Bout_ends(j)));
        T.BoutDisplacement(j) = sqrt(nansum(behaviorlog.f0_vx(T.Bout_starts(j):T.Bout_ends(j)))^2 +...
            nansum(behaviorlog.f0_vy(T.Bout_starts(j):T.Bout_ends(j)))^2);
        T.vigorSum(j) = sum(vigor(T.Bout_starts(j):T.Bout_ends(j)));
        T.vigorMean(j) = mean(vigor(T.Bout_starts(j):T.Bout_ends(j)));
        T.AngleMean(j) = mean(ang(T.Bout_starts(j):T.Bout_ends(j)));
        if ang(allpeaks(boutpeaks(1)))>0
            T.Laterality(j)=-1;
        else
            T.Laterality(j)=1;
        end
%         if j == 1
%             posBefore = 1;
%             [~, posAfter] = min(displint(T.Bout_ends(j):T.Bout_starts(j+1)));
%             posAfter = posAfter + T.Bout_ends(j);
%         elseif j == length(T.Bout_starts)
%             [~, posBefore] = min(displint(T.Bout_ends(j-1):T.Bout_starts(j)));
%             posBefore = posBefore + max(T.Bout_ends(j-1));
%             posAfter = length(displint);
%         else
%             [~, posBefore] = min(displint(T.Bout_ends(j-1):T.Bout_starts(j)));
%             posBefore = posBefore + max(T.Bout_ends(j-1));
%             [~, posAfter] = min(displint(T.Bout_ends(j):T.Bout_starts(j+1)));
%             posAfter = posAfter + T.Bout_ends(j);
%         end
        [x1, y1] = pol2cart(theta(T.Bout_starts(j)),1);
        [x2, y2] = pol2cart(theta(T.Bout_ends(j)),1);
        T.StartAngle(j) = theta(T.Bout_starts(j));
        T.BoutDirection(j) = atan2d(x1*y2-y1*x2,x1*x2+y1*y2);
        T.nankeep = nankeep;
        
        T.yaw(j) = max(abs(tmpYaw(T.Bout_starts(j):T.Bout_ends(j))));
        clear boutpeaks
    end
    TT(k).T = T;
    disp(num2str(k))
end
% %% plot pos color coded by speed
% % dp =[];
% % figure(1)
% % hold on
% % figure(2)
% % hold on
% % id1 = 1;
% % id2 = 1;
% for n = 3%length(TT)
%     figure
%     tiledlayout(2,2)
%     displacement = TT(n).T.displacement;
%     behaviorlog = rawData(n).behavior;
%     stimulus = rawData(n).stimulus;
%     tb = behaviorlog.t;
%     ts = stimulus.t;
%     vs = stimulus.conditional_velx;
%     vs(isnan(vs)) = 0;
%     [~,b,~] = unique(ts);
%     omrInt = interp1(ts(b), vs(b), tb, 'linear', 0);
%     omrOn = diff(omrInt)>0;
%     omrOff = diff(omrInt)<0;
%     for k = [0 3 6 9]
%         if k == 0 
%             stimstart = find(omrOff,1);
%             stimend = find(omrOn);
%             stimend = stimend(find(stimend>stimstart,1));
%         else
%         stimstart = find(omrInt==k,1);
%         stimend = find(omrOff);
%         stimend = stimend(find(stimend>stimstart,1));
%         end
%     
%     cm = cmocean('speed', 100);
%     displint = interp1(find(~isnan(displacement)), displacement(~isnan(displacement)), 1:length(displacement), 'linear', 0);
%     displint = medfilt1(displint,7);
%     
% %     scatter(behaviorlog.f0_x, behaviorlog.f0_y,nrm(displint).*30+0.1,cm(round(nrm(displint).*99)'+1,:), 'filled')
% % scatter(smooth(behaviorlog.f0_x, 500, 'moving'), smooth(behaviorlog.f0_y, 500, 'moving'),5,cm(round(nrm(displint).*99)'+1,:), 'filled')
% nexttile
% % scatter(pxscale.*behaviorlog.f0_x(stimstart:stimend), pxscale.*behaviorlog.f0_y(stimstart:stimend),'k', 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none')
% scatter(pxscale.*behaviorlog.f0_x(omrInt==k & behaviorlog.t>600), pxscale.*behaviorlog.f0_y(omrInt==k & behaviorlog.t>600),'k', 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none')
% 
% % plot(behaviorlog.f0_x, behaviorlog.f0_y, 'Color', [0 0 0 .3])
% % daspect([1 1 1])
% % x = behaviorlog.f0_x;
% % y = behaviorlog.f0_y;
% % x = smooth(x, 1, 'moving');
% % y = smooth(y, 1, 'moving');
% % x(movstd(x,25)>30) = NaN;
% % y(movstd(y,25)>30) = NaN;
% % figure
% % plot(x)
% % hold on
% % plot(y)
% title([master.knockout(n)])
% x = smooth(behaviorlog.f0_x, 1, 'moving');
% y = smooth(behaviorlog.f0_y, 1, 'moving');
% % if master.knockout(n)=='+'
% %     figure(1)
% %     subplot(4,4,id1)
% % % scatter(x, y,5,'r', 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none')
% % p= plot(x,y, 'k', 'LineWidth', .5);
% % % p.Color(4) = 0.2;
% % id1 = id1 +1;
% % else 
% %     figure(2)
% %     subplot(4,4,id2)
% % %     scatter(x, y,5,'r', 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none')
% % p= plot(x,y,  'k', 'LineWidth', .5);
% % % p.Color(4) = 0.2;
% % id2 = id2 +1;
% % end
% dp(n) = nansum(sqrt(diff(x).^2 + diff(y).^2));
%     % plot(behaviorlog.f0_x, behaviorlog.f0_y,'k')
% %     hold on
% %     % scatter(behaviorlog.f0_x, behaviorlog.f0_y,4,cm(round(nrm(displint).*99)'+1,:), 'filled')
%     daspect([1 1 1])
%     xlim([0 100])
%     ylim([0 50])
% % %     colormap(cm)
% % %     set(gca, 'Color', [0.7 0.7 0.7])
% % %     set(gcf, 'Color', [0.7 0.7 0.7])
% %     % set(gca, 'Visible', 'off')
% %     
%     xticks([0 700])
%     yticks([0 150])
% %     title([master.knockout(n) num2str(dp(n))])
% % title(['#' num2str(n)])
%     end
% end
% % figure
% % distIdx = master.knockout=='+';
% % plotSpread(dp, 'distributionIdx', distIdx)
%%
rat = [];
for n = 1:length(TT)
%     figure
%     plot(TT(n).T.ang)
%     hold on
%     plot(diff(TT(n).T.displacement))
%     plot(TT(n).T.displacement)
%     plot(TT(n).T.theta)
% plot(medfilt1(diff(TT(n).T.theta),3))
%     scatter(TT(n).T.allPeaks, TT(n).T.ang(TT(n).T.allPeaks), 'k', 'filled', 'SizeData', 3)
%     scatter(TT(n).T.Bout_starts,TT(n).T.ang(TT(n).T.Bout_starts), 'g', 'filled', 'SizeData', 3) 
%     scatter(TT(n).T.Bout_ends,TT(n).T.ang(TT(n).T.Bout_ends), 'r', 'filled', 'SizeData', 3)
    sc = zeros(length(TT(n).T.Bout_starts), 1);
    for bs = 1:length(TT(n).T.Bout_starts)
        sc(bs) = std(1./diff(TT(n).T.allPeaks(TT(n).T.allPeaks>TT(n).T.Bout_starts(bs) & TT(n).T.allPeaks<TT(n).T.Bout_ends(bs))));
        %         if sc(bs)>0.07
        %             plot(TT(n).T.Bout_starts(bs):TT(n).T.Bout_ends(bs), TT(n).T.ang(TT(n).T.Bout_starts(bs):TT(n).T.Bout_ends(bs)), 'r')
        %         end
        if TT(n).T.Max_amp(bs) > .7 && TT(n).T.BoutDirection(bs)>-10 && TT(n).T.BoutDirection(bs)<10
%             plot(TT(n).T.Bout_starts(bs):TT(n).T.Bout_ends(bs), TT(n).T.ang(TT(n).T.Bout_starts(bs):TT(n).T.Bout_ends(bs)), 'g')
        elseif TT(n).T.Max_amp(bs) < .7 && TT(n).T.BoutDirection(bs)>-10 && TT(n).T.BoutDirection(bs)<10
%             plot(TT(n).T.Bout_starts(bs):TT(n).T.Bout_ends(bs), TT(n).T.ang(TT(n).T.Bout_starts(bs):TT(n).T.Bout_ends(bs)), 'r')
        end
    end
    TT(n).T.freqStd = sc;
%     scatter(TT(n).T.Bout_starts, sc.*10, 'm', 'filled')
    rat = cat(1,  rat, [length(TT(n).T.Bout_starts) sum(sc>0.07)]);
%     plot(TT(n).T.allPeaks, [0 diff(TT(n).T.allPeaks)], 'm')
% title(master.knockout(n))
end
dontKeep = find((rat(:,2)./rat(:,1))>0.1);
% %% plot pos color coded by fish speed and dot size is omr speed
% cm = cmocean('speed', 100);
% displint = interp1(find(~isnan(displacement)), displacement(~isnan(displacement)), 1:length(displacement));
% displint = medfilt1(displint,7);
% tb = behaviorlog.t;
% ts = stimuluslog.t;
% vs = stimuluslog.conditional_velx;
% vs(isnan(vs)) = 0;
% [a,b,c] = unique(ts);
% omrInt = interp1(ts(b), vs(b), tb, 'linear', 0);
% figure
% scatter(behaviorlog.f0_x, behaviorlog.f0_y,nrm(omrInt).*30+0.1,cm(round(nrm(displint).*99)'+1,:), 'filled')
% % plot(behaviorlog.f0_x, behaviorlog.f0_y,'k')
% hold on
% % scatter(behaviorlog.f0_x, behaviorlog.f0_y,4,cm(round(nrm(displint).*99)'+1,:), 'filled')
% daspect([1 1 1])
% xlim([0 1504])
% colormap(cm)
% set(gca, 'Color', [0.7 0.7 0.7])
% set(gcf, 'Color', [0.7 0.7 0.7])
% % set(gca, 'Visible', 'off')
% xticks([0 1504])
% yticks([0 350])
% %% plot tail angle and bout start/end
% for n = 1
%     behaviorlog = rawData(n).behavior;
%     T = TT(n).T;
%     plotAng = TT(n).T.ang;
%     stimulus = rawData(n).stimulus;
%     figure
%     hold on
%     plot(behaviorlog.t, plotAng, 'k')
%     plot(stimulus.t, stimulus.conditional_velx,'b')
%     
%     pks = vertcat(T.Boutpeaks{:});
%     scatter(behaviorlog.t(pks), plotAng(pks), 'k', 'filled')
%     scatter(behaviorlog.t(T.Bout_starts), plotAng(T.Bout_starts), 'g', 'filled')
%     scatter(behaviorlog.t(T.Bout_ends), plotAng(T.Bout_ends), 'r', 'filled')
%     % plot(behaviorlog.t, omrInt./10, 'r')
%     % scatter(behaviorlog.t(T.Bout_starts), T.BoutDirection./100, 'm', 'filled')
%     % scatter(behaviorlog.t(T.Bout_starts), T.AngleMean, 'y', 'filled')
% end
%% make a long data table with each bout

sw = 0;
for n = 1:length(TT)
    sw = sw+length(TT(n).T.Bout_starts);
end
longTable = table('Size', [sw, 30], 'variableTypes', {'double', 'double',...
    'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double',...
    'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'categorical'}, 'variableNames', {'Duration', 'Ncycles', 'TBF',...
    'MaxAmp', 'MeanAmp', 'MeanSpeed', 'VigorSum', 'VigorMean', 'TotalVigor', 'BoutSpeed', 'TotalSpeed', 'BoutDirection', 'BoutDisplacement', 'AngleMean', 'Yaw',...
    'dCoeffVar', 'CoeffVar', 'BoutPeaks', 'Distance', 'fishID', 'StimVel', 'InterboutDuration', 'BoutId', 'FreqStd', 'PeriodStd', 'dPeriodStd', 'StartAngle', 'omrDirection', 'XStart', 'Genotype'});
id = 1;
for n = 1:length(TT)
    for k = 1:length(TT(n).T.Bout_starts)
        longTable.Duration(id) = TT(n).T.Bout_duration(k)/1000;
        longTable.Ncycles(id) = TT(n).T.N_cycles(k);
        longTable.TBF(id) = TT(n).T.TBF(k);
        longTable.MaxAmp(id) = TT(n).T.Max_amp(k);
        longTable.MeanSpeed(id) = TT(n).T.MeanSpeed(k);
        longTable.BoutDirection(id) = TT(n).T.BoutDirection(k);
        longTable.MeanAmp(id) = TT(n).T.Mean_amp(k);
        longTable.AngleMean(id) = TT(n).T.AngleMean(k);
        longTable.fishID(id) = n;
        longTable.StimVel(id) = TT(n).T.StimVel(k);
        longTable.Genotype(id) = master.knockout(n);
        longTable.Distance(id) = TT(n).T.distance(k);
        longTable.StartAngle(id) = TT(n).T.StartAngle(k);
        longTable.XStart(id) = TT(n).T.Xpos_start(k);
        longTable.BoutDisplacement(id) = TT(n).T.BoutDisplacement(k);
%         longTable.protN(id) = TT(n).T.prodN(k);
        longTable.BoutId(id) = k;
        longTable.omrDirection(id) = TT(n).T.omrDirection(k);
        longTable.FreqStd(id) = TT(n).T.freqStd(k);
        longTable.PeriodStd(id) = TT(n).T.periodStd(k);
        longTable.dPeriodStd(id) = TT(n).T.dperiodStd(k);
        longTable.CoeffVar(id) = longTable.PeriodStd(id)/TT(n).T.meanPeriod(k);
        longTable.dCoeffVar(id) = longTable.dPeriodStd(id)/TT(n).T.meandPeriod(k);
        longTable.Yaw(id) = TT(n).T.yaw(k);
        longTable.InterboutDuration(id) = TT(n).T.Interbout_duration(k)/1000;
        longTable.BoutSpeed(id) = longTable.BoutDisplacement(id)/longTable.Duration(id);
        longTable.TotalSpeed(id) = longTable.BoutDisplacement(id)/(longTable.Duration(id) + longTable.InterboutDuration(id));
        longTable.VigorSum(id) = TT(n).T.vigorSum(k);
        longTable.VigorMean(id) = TT(n).T.vigorMean(k);
        longTable.TotalVigor(id) = longTable.VigorSum(id)/(longTable.Duration(id) + longTable.InterboutDuration(id));
        id = id + 1;
        if mod(id, 100)==0
            disp(num2str(id))
        end
    end
end

% %% total number of bouts/fish
% scd = [];
% sci = [];
% genotp = [];
% 
% for n = 1:length(TT)
%     if ~ismember(n, dontKeep)
%         disp(num2str(n))
%         %     displacement = TT(n).T.displacement;
%         %     displint = interp1(find(~isnan(displacement)), displacement(~isnan(displacement)), 1:length(displacement), 'linear', 0);
%         %     displint = medfilt1(displint,7);
%         %     scd = [scd sum(displint)./length(displint)]
%         % scd = [scd length(TT(n).T.Bout_starts)];
%         tmpsci = zeros(4,1);
%         tmpgenotp = zeros(4,1);
%         tmpscd = zeros(4,1);
%         els = zeros(4,1);
%         for k = 0:3:9
%             %
% %             selector =  longTable.Duration<400 & longTable.fishID == n & longTable.StimVel == k &...
% %                 longTable.Ncycles<30 &  longTable.AngleMean < 0.2 &...
% %                 longTable.TBF>15 & longTable.MeanAmp > 0.3 ;
%             %
%             
%                                 selector =  longTable.fishID == n & longTable.StimVel == k & bsLong>300*60*10; %ignore first 10 minutes
%          
%             
%             %    scd = [scd sum(longTable.Distance(selector))];
%             behaviorlog = rawData(n).behavior;
%             stimuluslog = rawData(n).stimulus;
%             tb = behaviorlog.t;
%             ts = stimuluslog.t;
%             vs = stimuluslog.conditional_velx;
%             vs(isnan(vs)) = 0;
%             [~,b,~] = unique(ts);
%             omrInt = interp1(ts(b), vs(b), tb, 'linear', 0);
%             ignore = false(size(omrInt));
%             ignore(300*10*60:end) = true;
%             stimTime = sum(abs(omrInt(TT(n).T.nankeep==0 & ignore))==k)/framerate;
%             tmpscd(k/3 +1) = length(longTable.Distance(selector))./stimTime;
%             els(k/3 +1) = length(longTable.Distance(selector));
%             
%             tmpsci(k/3 +1) = k;
%             tmpgenotp(k/3 +1) = master.knockout(n) == '+';
%         end
%         if sum(els>20)==4
%             scd = [scd tmpscd'];
%             sci = [sci tmpsci'];
%             genotp = [genotp tmpgenotp'];
%         end
%     end
% end
% sci(genotp==1) = sci(genotp==1)+1;
% figure
% hold on
% plotSpread(scd(genotp==0), 'distributionId', sci(genotp==0), 'distributionColor', 'k')
% hold on
% plotSpread(scd(genotp==1), 'distributionId', sci(genotp==1), 'distributionColor', 'r')
% sci(genotp==1) = sci(genotp==1)-1;
% % plotSpread(scd, 'distributionIdx', sci(ismember(sci, 0:3:18)), 'xValues', 0:3:18, 'distributionColors', 'r' )
% % plotSpread(scd, 'distributionIdx', sci(ismember(sci, 21:3:39)), 'xValues', 1:3:19, 'distributionColors', 'k')
% % hold on
% % boxchart(sci(ismember(sci, 0:3:18)), scd(ismember(sci, 0:3:18)), 'BoxFaceColor', 'r', 'MarkerColor', 'r', 'MarkerStyle', '.')
% % boxchart(sci(ismember(sci, 21:3:39))-20, scd(ismember(sci, 21:3:39)), 'BoxFaceColor', 'k', 'MarkerColor', 'k', 'MarkerStyle', '.')
% val = zeros(4,1);
% for n = 1:4
%     val(n) = mean(scd(genotp==0 & sci==(n-1)*3));
%     sem(n) = std(scd(genotp==0 & sci==(n-1)*3))/sqrt(numel(scd(genotp==0 & sci==(n-1)*3)));
% end
% errorbar(0:3:9, val,sem, 'k')
% val = zeros(4,1);
% for n = 1:4
%     val(n) = mean(scd(genotp==1 & sci==(n-1)*3));
%     sem(n) = std(scd(genotp==1 & sci==(n-1)*3))/sqrt(numel(scd(genotp==1 & sci==(n-1)*3)));
% end
% errorbar(0:3:9, val,sem, 'r')
% 
% xlabel('OMR speed')
% ylabel('bout rate')
% ft = table(scd(:), sci(:), genotp(:), 'VariableNames', {'boutRate', 'omr', 'genotype'});
% lm1 = fitlm(ft, 'boutRate ~ omr*genotype')

%% calculate time per fish spend per omr speed
totTime = zeros(length(TT), length(0:3:18));
speeds = 0:3:18;
framerate = 300;
for n = 1:length(TT)
    badSwims = zeros(size(TT(n).T.ang));
    for sws = 1:length(TT(n).T.Bout_starts)
        if TT(n).T.freqStd > 0.07
            badSwims(TT(n).T.Bout_starts(sws):TT(n).T.Bout_ends(sws)) = 1;
        end
    end
            
    for omr = 1:7
        behaviorlog = rawData(n).behavior;
        stimuluslog = rawData(n).stimulus;
        tb = behaviorlog.t;
        ts = stimuluslog.t;
        vs = stimuluslog.conditional_velx;
        vs(isnan(vs)) = 0;
        [~,b,~] = unique(ts);
        omrInt = interp1(ts(b), vs(b), tb, 'linear', 0);
        stimTime = sum(abs(omrInt(TT(n).T.nankeep==0 & badSwims'==0))==speeds(omr))/framerate;
        totTime(n,omr) = stimTime;
    end
end


%% fish average per omr speed
figure
tiledlayout(2,3,'Padding', 'tight', 'TileSpacing', 'tigh')
metric = {'Duration', 'InterboutDuration', 'BoutSpeed', 'MeanAmp', 'TBF', 'TotalSpeed'};
for met = 1:length(metric)
    nexttile
    hold on
    mKO = [];
    mWT = [];
    
        val = [];
        els = [];
        distId = [];
        idx = 1;
        speedSel = {0, 3, 6, 9};
    for omr = 1:4
        sel = unique(longTable.fishID);
        sel(ismember(sel, dontKeep)) = [];
%         sel(sel==18) = [];
        for n = 1:length(sel)  %& (longTable.BoutDirection >-20 & longTable.BoutDirection <20) & omrFollow
            if omr > 1
                val(idx,n) = nanmedian(eval(['longTable.' metric{met} '(longTable.fishID==sel(n)  & longTable.FreqStd <= 0.07   & ~ismember(longTable.fishID, dontKeep)  & ismember(longTable.StimVel, speedSel{' num2str(omr) '}))']));
                els(idx,n) = length(eval(['longTable.' metric{met} '(longTable.fishID==sel(n)  & longTable.FreqStd <= 0.07   & ~ismember(longTable.fishID, dontKeep)  & ismember(longTable.StimVel, speedSel{' num2str(omr) '}))']));
%                 val(idx,n) = val(idx,n)./sum(totTime(sel(n),1 + speedSel{omr}/3));
%                 val(idx,n) = val(idx,n)./(sum(posTrial(sel(n), 1 + speedSel{omr}/3))/framerate);
%                 val(idx,n) = val(idx,n)./val(1,n);
            else
                val(idx,n) = nanmedian(eval(['longTable.' metric{met} '(longTable.fishID==sel(n)  & longTable.FreqStd <= 0.07   & ~ismember(longTable.fishID, dontKeep) & ismember(longTable.StimVel, speedSel{' num2str(omr) '}))']));
                els(idx,n) = length(eval(['longTable.' metric{met} '(longTable.fishID==sel(n)  & longTable.FreqStd <= 0.07   & ~ismember(longTable.fishID, dontKeep) & ismember(longTable.StimVel, speedSel{' num2str(omr) '}))']));
%                 val(idx,n) = val(idx,n)./sum(totTime(sel(n),1 + speedSel{omr}/3));
%                 val(idx,n) = val(idx,n)./sum(posTrial(sel(n), 1 + speedSel{omr}/3));
%                 val(idx,n) = val(idx,n)./((300*600+300*5*10*6)/framerate);
            end
            if master.knockout(sel(n)) == '-'
                distId(n) = 0;
            else
                distId(n) = 1;
            end
        end
        idx = idx + 1;
    end
elss = els>20;
if strcmp(metric{met}, 'BoutSpeed') || strcmp(metric{met}, 'TotalSpeed')
    val = val.*pxscale;
end
% val = (val-nanmean(val))./nanmean(val);
errorbar(0:3:9, mean(val(:,distId==0 & sum(elss,1)==4),2),std(val(:,distId==0 & sum(elss,1)==4),[],2)./sqrt(size(val(:,distId==0 & sum(elss,1)==4),2)), 'k')
errorbar(0:3:9, mean(val(:,distId==1 & sum(elss,1)==4),2),std(val(:,distId==1 & sum(elss,1)==4),[],2)./sqrt(size(val(:,distId==1 & sum(elss,1)==4),2)), 'r')
ddistId = repmat(distId(sum(elss,1)==4), 4, 1);
ddistId = ddistId(:);
xv = repmat([0:3:9]', 1, size(val(:,sum(elss,1)==4),2));
xv(ddistId==1) = xv(ddistId==1)+1;
xv = xv(:);

disp(['wt: ' num2str(mean(mean(val(:,distId==0 & sum(elss,1)==4),1)))...
    ' ± ' num2str(std(mean(val(:,distId==0 & sum(elss,1)==4),1))/sqrt(sum(distId==0)))])

disp(['ko: ' num2str(mean(mean(val(:,distId==1 & sum(elss,1)==4),1)))...
    ' ± ' num2str(std(mean(val(:,distId==1 & sum(elss,1)==4),1))/sqrt(sum(distId==1)))])

val2 = val(:,sum(elss,1)==4);
fid = zeros(size(val2));
fid = repmat([1:size(fid,2)],size(fid,1),1);
fid = fid(:);
val2 = val2(:);
plotSpread(val2(ddistId==0), 'distributionId', xv(ddistId==0), 'distributionColor', 'k')
hold on
plotSpread(val2(ddistId==1), 'distributionId', xv(ddistId==1), 'distributionColor', 'r')



title(metric{met})
currylim = get(gca, 'Ylim');
ylim([0 currylim(2)])
end
ft = table(val2, xv, ddistId, fid, 'VariableNames', {'distance', 'omr', 'genotype', 'fishId'});
lm1 = fitlme(ft, 'distance ~ omr*genotype + (omr|fishId)');
%% overall bout rate
 sel = unique(longTable.fishID);
 sel(ismember(sel, dontKeep)) = [];
 sel(sum(elss,1)~=4) = [];
 nb = zeros(size(sel));
 totDur = zeros(size(sel));
 gt = zeros(size(sel));
 for n = 1:length(sel)
%      nb(n) = length(longTable.Duration(longTable.fishID==sel(n)  & longTable.FreqStd <= 0.07 & ismember(longTable.StimVel,[0 3 6 9])));
     nb(n) = nanmedian(longTable.CoeffVar(longTable.fishID==sel(n)  & longTable.FreqStd <= 0.07 & ismember(longTable.StimVel,[0 3 6 9])));
%      nb(n) = nanmedian(longTable.FreqStd(longTable.fishID==sel(n)  & longTable.FreqStd <= 0.07 & ismember(longTable.StimVel,[0 3 6 9]))./...
%          longTable.TBF(longTable.fishID==sel(n)  & longTable.FreqStd <= 0.07 & ismember(longTable.StimVel,[0 3 6 9])));
     totDur(n) = height(rawData(sel(n)).behavior)/framerate;
     
     if master.knockout(sel(n)) == '-'
         gt(n) = 1;
     else
         gt(n) = 2;
     end
 end
%  br = nb./totDur;
br = nb;
figure
hold on
plotSpread(br, 'distributionId', gt, 'distributionColor', 'k')
av = [nanmean(br(gt==1)) nanmean(br(gt==2))]
sem = [std(br(gt==1))/sqrt(sum(gt==1)) std(br(gt==2))/sqrt(sum(gt==2))]
errorbar(unique(gt), av, sem, 'r')
yl = get(gca, 'Ylim');
ylim([0 yl(2)])
% ylabel('Bout rate (Hz)')
ylabel('Coefficient of variation half-cycle period')
[h,p] = ttest2(br(gt==1), br(gt==2))
text(1.2,0.2,['p = ' num2str(p)])
xticklabels({'WT', 'KO'})

 


% %% get bout starts in long table format
% bsLong = zeros(height(longTable),1);
% idx = 1;
% for n = 1:length(TT)
%     for k = 1:length(TT(n).T.Bout_starts)
%         bsLong(idx) = TT(n).T.Bout_starts(k);
%         idx = idx + 1;
%     end
% end


