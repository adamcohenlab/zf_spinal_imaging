%iterate over all ica_pca_selected files in folder
files = dir('*ica_pca_selected.mat');
for fileIter = 1:length(files)
    fileN = str2num(cell2mat(regexp(files(fileIter).name,'\d*','match')));
    lN = str2num(cell2mat(regexp(pwd,'(?<=l)\d*','match')));
    if exist(['l' num2str(lN) '-tr' num2str(fileN) '_VNR.mat'], 'file')
        % load VNR file
        load(['l' num2str(lN) '-tr' num2str(fileN) '_VNR.mat']);
        % load mean frame
        imgMat = loadtiff(strrep(files(fileIter).name, 'ica_pca_selected.mat', 'meanFrame.tif'));
        % load rois
        load(strrep(files(fileIter).name, 'ica_pca_selected', 'rois'));
        % load traces
        load(files(fileIter).name)
        
        % hack to fill empty traces, will only work if the first one is not
        % empty
        for n = 1:length(trace)
            if isempty(trace{n})
                trace{n} = zeros(size(trace{1}));
            end
        end
        %get traces
        tr = cell2mat(trace);
        % traces are chopped by 500 frames, data in burstImage is not
        ba = burstImg-499;
        % bursts and instantaneous freq outside of imaging time are removed
        noBurst = ba<1 | ba>size(tr,1);
        ba(noBurst) = [];
        instFreq(noBurst) = [];
        % initialize variables
        phb = [];
        idx = 1;
        indPeak = [];
        indVal = [];
        indTangVal = [];
        indTangPeak = [];
        %loop over bursts
         for n = 1:length(ba)-1
             %only use bursts that lie within 15 to 40 inst freq range
            if instFreq(n+1)>15 && instFreq(n+1)<40
                % frames of the current cycle (burst to burst)
                interval = ba(n):ba(n+1);
                newInterval = linspace(interval(1),interval(end),100);
                % interpolate trace to havce 100 saples
                tmpTr = interp1(interval,tr(interval,:),newInterval,'pchip');
                % loop over cells
                for k = 1:size(tmpTr,2)
%                     tmpTr(:,k) = smooth(tmpTr(:,k),10,'moving');
                    % get index of of the maximum value
                    [indVal(idx,k),indPeak(idx,k)] = max(tmpTr(:,k));
                end
                
                phb(:,:,idx)  = tmpTr;

                idx = idx + 1;
            end
        end
  
        
        av = squeeze(mean(phb,3));
        a = rois;
        clear sbar
        cc = cmocean('phase',100);
        sbar(:,:,1) = repmat(cc(:,1),1,10)';
        sbar(:,:,2) = repmat(cc(:,2),1,10)';
        sbar(:,:,3) = repmat(cc(:,3),1,10)';
        timf = zeros(size(imgMat,1), size(imgMat,2), 3);
        idx = 1;
        idx2 = 0;
        gain = 1;
        [v,~] = max(av);
        % v = mean(indVal,1);
        % convert to complex number and average to get average angle
        pk2 = angle(mean(exp(1i*pi*(indPeak./50)),1))./pi;
%         pkTang = angle(mean(exp(1i*pi*(indTangPeak./50)),1))./pi;
        pk2 = mod((pk2+2),2)*50;
%         pkTang = mod((pk2+2),2)*50;
        y = pk2;
%         y = pkTang;
        y(pk2<25) = y(pk2<25)+50;
        y(pk2>=75) = y(pk2>=75)-50;
        pk = ceil(pk2);
        
        interpStdV = interp1(1/50000:1/50000:25,stdV,1/996.4:1/996.4:25);
         smTr = [];
        xc = [];
        for n = 1:size(tr,2)
            smTr(:,n) = smooth(tr(:,n),10,'moving');
        end
        for n = 1:size(tr,2)
            xc(:,n) = xcorr(interpStdV(500:end),smTr(:,n),50);
        end
        [~,pkXc] = max(xc);
        
       
        [pw,px] = pspectrum(interpStdV,996.4);
%         pw = mean(pw,2);
        [~,fIdx] = max(pw(px>15));
        f  = px(px>15);
        f = f(fIdx);
        
        % av(:,pk2<25 | pk2>=75) = circshift(av(:,pk2<25 | pk2>=75),50,1);
        % pk(pk<25) = pk(pk<25)+50;
        % pk(pk>=75) = pk(pk>=75)-50;
        % pk = pk-25;
        % cmm = hsv(size(av,1));
        % cm = cmm(p,:);
        % % cm = lines(size(tr,2));
        % cm = brighten(cm,0);
        % cm = flipud(cm);
        % cmm = cmocean('phase',size(av,1));
        cmm = cmocean('phase',100);
        % cmm = hsv(size(av,1));
        cm = cmm(pk,:);
        value = mean(imgMat,3);
        value = value - min(value(:));
        value = value./(max(value(:))*0.5);
        value = ones(size(value));
        normalization = zeros(size(value));
        xsort = [];
        for n = 1:length(space)
            tmp = space{n};
            for p = 1:size(tmp,2)
                if ~ismember(idx,del)
                    idx2 = idx2 + 1;
                    %             if idx2 == 8 || idx2 == 9
                    %             ymin = round(min(a{1,n}(:,1)));
                    %             ymax = round(max(a{1,n}(:,1)));
                    %             xmin = round(min(a{1,n}(:,2)));
                    %             xmax = round(max(a{1,n}(:,2)));
                    ymin = round(min(a{1,n}(:,1)));
                    ymax = min([round(max(a{1,n}(:,1))) size(imgMat,2)]);
                    xmin = round(min(a{1,n}(:,2)));
                    xmax = min([round(max(a{1,n}(:,2))) size(imgMat,1)]);
                    xsort(n) = ymin;
                    ic = tmp(:,p);
                    ic = ic-min(ic);
                    ic = ic./max(ic);
                    icm = flipud(toimg(ic, (xmax-xmin)+1, (ymax-ymin)+1));
                    saturation = zeros(size(imgMat,1), size(imgMat,2));
                    saturation(xmin:xmax,ymin:ymax) = icm;
                    if v(idx2)<0.5
                        cm(idx,:) = 1;
                    end
                    valuetmp = value;
                    valuetmp(saturation(:)==0) = 0;
                    tim = zeros(size(value,1), size(value,2), 3);
                    hh = rgb2hsv(cm(idx2,:));
                    %             tim(:,:,1) = hh(1);
                    %             tim(:,:,2) = valuetmp.*gain;
                    %             tim(:,:,3) = saturation;
                    tim(:,:,1) = saturation.*cm(idx2,1);
                    tim(:,:,2) = saturation.*cm(idx2,2);
                    tim(:,:,3) = saturation.*cm(idx2,3);
                    normalization = normalization + logical(valuetmp);
                    timf = timf + tim;
                    %             end
                end
                idx = idx + 1;
            end
        end
        timf(end-9:end,end-99:end,:) = sbar;
        figure(444); clf
        subplot(2,1,1)
        imshow2(timf)
        subplot(2,1,2)
        set(gca,'Color','k')
        set(gcf,'Color','k')
        hold on
        allTr = cell2mat(trace);
        tmpPk = pk;
        tmpPk(v<0.5) = 0;
        % figure
        % hold on
        [val,idx] = sort(tmpPk);
        % [~,idx] = sort(xsort);
        for n = 1:size(allTr,2)
            plot(allTr(:,idx(n))+n*1.3*(max(allTr(:)) - min(allTr(:))),'Color', cm(idx(n),:))
            %     plot(allTr(:,idx(n))+val(n)*5+n*4,'Color', cm(idx(n),:))
        end
        saveas(gcf, strrep(files(fileIter).name, 'ica_pca_selected.mat', 'VNR_trig_phase.fig'));
        save(strrep(files(fileIter).name, 'ica_pca_selected.mat', 'VNR_trig_phase.mat'),'av', 'phb', 'indPeak', 'pk2', 'xsort',...
            'xc', 'pkXc', 'f');
        
        figure(555); clf
        % hold on
        xs = xsort/1024;
        xs = xs*100;
        xs = round(xs);
        xs(xs==0) = 1;
        for n = 1:length(xs)
            polarplot(linspace(0,2*pi,100),av(:,n),'Color',cmm(xs(n),:))
            hold on
        end
        saveas(gcf, strrep(files(fileIter).name, 'ica_pca_selected.mat', 'VNR_trig_polar.fig'));
    end
end
%%
[~,avTang] = max(diff(av));
x = xsort;
% y = pk2;
% y(pk2<25) = y(pk2<25)+50;
% y(pk2>=75) = y(pk2>=75)-50;
% y = pkTang;
% y(pkTang<25) = y(pkTang<25)+50;
% y(pkTang>=75) = y(pkTang>=75)-50;
% y = avTang;
% y(avTang<25) = y(avTang<25)+50;
% y(avTang>=75) = y(avTang>=75)-50;
y = pkXc-50;
yy = pkXc-50;
y(yy<-12.5) = y(yy<-12.5)+17.5;
y(yy>=5) = y(yy>=5)-17.5;
figure
scatter(x(max(av)>.5),y(max(av)>.5),'filled')
%%
[~,idx] = sort(xsort);
figure
stackplot(tr(:,idx))
%%
smTr = [];
xc = [];
for n = 1:size(tr,2)
    smTr(:,n) = smooth(tr(:,n),10,'moving');
end
ref = idx(max(av)>.5);
for n = 1:size(tr,2)
%     xc(:,n) = xcorr(tt(500:end),smTr(:,n),50);
xc(:,n) = xcorr(smTr(:,ref(1)),smTr(:,n),50);
end
figure
stackplot(smTr)
figure
plot(xc)
%%
% [~,pkXc] = max(xc);
% y(avTang<25) = y(avTang<25)+50;
% y(avTang>=75) = y(avTang>=75)-50;
% y = pkXc-50;
y = px(pkXc)-px(50);
yy = y;
y(yy<-12.5) = y(yy<-12.5)+(996.4/f)/2;
y(yy>=5) = y(yy>=5)-(996.4/f)/2;
figure
hold on
for n = 1:length(y)
    if cy(n) > (cx(n)*fobj.p1 + fobj.p2) 
scatter(cx(n),y(n),'r', 'filled')
    else
        scatter(cx(n),y(n),'b', 'filled')
    end
end