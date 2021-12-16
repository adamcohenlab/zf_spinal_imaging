%%load data

data = [];
stdVAll = [];
% master = omrspeedgainmaster;

xyPrev = 0;
zPrev = 0;
idx1 = 1;
for fileIter = 1:height(master)
    folder = master.folder{fileIter};  
    file = master.fileName{fileIter};
    path = [folder '\' file];
    fileN = str2num(cell2mat(regexp(file,'\d*','match')));
    lN = str2num(cell2mat(regexp(folder,'(?<=\\l)\d*','match')));
    load(path)
    load(strrep(path,'ica_pca_selected', 'rois'))
    load([folder '\l' num2str(lN) '-tr' num2str(fileN) '_VNR.mat']);
    load(strrep(path,'ica_pca_selected', 'VNR_trig_phase2'))

     if xyPrev ~= master.xyposId(fileIter) || zPrev ~= master.zPos(fileIter)
        trial = 1;
    else
        trial = trial + 1;
     end
     stdVAll{fileIter} = stdV;
    for iter2 = 1:size(tr,2)
        data(idx1).trace = tr(:,iter2)./std(tr(:,iter2));
        data(idx1).burstImg = burstImg;
        data(idx1).boutEndIdx = boutEndIdx;
        data(idx1).boutStartIdx = boutStartIdx;
        data(idx1).rois = rois(iter2);
        data(idx1).instFreq = instFreq;
        data(idx1).stdV = fileIter;
        data(idx1).burstAmp = stdV(bursts);
        data(idx1).fishId = master.fishId(fileIter);
        data(idx1).xypos = master.xyposId(fileIter);
        data(idx1).brainRegion = master.brainRegion(fileIter);
        data(idx1).gainAdjust = master.gainAdjust(fileIter);
        data(idx1).omrGain = master.omrGain(fileIter);
        data(idx1).omrSpeed = master.omrSpeed(fileIter);
        data(idx1).zPos = master.zPos(fileIter);
        data(idx1).trial = trial;
        xyPrev = master.xyposId(fileIter);
        zPrev = master.zPos(fileIter);
        idx1 = idx1 + 1;
    end
     
end
save('complete_dataset.mat', 'data')
save('complete_stdVAll.mat', 'stdVAll')
%%
clear all
load('complete_dataset.mat')
%%
for n = 1:length(data)
    if sum(isnan(data(n).trace))==0
    [av, pk2, pkXc, f, indivInstFreq, indCycle] = calcVNRphase(data(n).burstImg, stdVAll{data(n).stdV},...
        data(n).trace, data(n).instFreq);
    else
        av = [];
        pkXc = [];
        f = [];
        pk2 = [];
        indivInstFreq = [];
        indCycle = [];
    end
    data(n).av = av;
    data(n).pkXc = pkXc;
    data(n).f = f;
    data(n).pk2 = pk2;
    data(n).selInstFreq = indivInstFreq;
    data(n).indCycle = indCycle;
    if mod(n,10)==0
        disp(['cell ' num2str(n) ' of ' num2str(length(data))])
    end
end
save('complete_dataset.mat', 'data')
save('complete_stdVAll.mat', 'stdVAll')
