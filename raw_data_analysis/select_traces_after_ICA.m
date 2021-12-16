trace_files = dir('*trace*.mat');
space_files = dir('*space*.mat');

for k = 1:length(trace_files)
    if ~exist(strrep(trace_files(k).name, 'trace', 'selected'), 'file')
        
        load(trace_files(k).name);
        load(space_files(k).name);
        load(strrep(trace_files(k).name, 'ica_pca_trace.mat', 'rois.mat'));
        a = rois;
        img = loadtiff(strrep(trace_files(k).name, 'ica_pca_trace.mat', 'meanFrame.tif'));
%         load(strrep(trace_files(k).name, 'ica_pca_trace.mat', 'meanFrame.tif'));
        for n =1:size(trace,2)
            ymin = round(min(a{1,n}(:,1)));
            ymin = max([ymin 1]);
            ymax = min([round(max(a{1,n}(:,1))) size(img,2)]);
            xmin = round(min(a{1,n}(:,2)));
            xmin = max([xmin 1]);
            xmax = min([round(max(a{1,n}(:,2))) size(img,1)]);
            M1 = img(xmin:xmax,ymin:ymax);
            icsTime = trace{n};
            maxFrame = size(icsTime,1);
            trTmp = [];
            trDisp = [];
            for m = 1:size(icsTime,2)
                %                 tmp = icsTime(:,m) - medfilt1(icsTime(:,m),2001);
                try
                    tmp = correct_bleaching(icsTime(200:maxFrame-round(maxFrame/100),m));
                catch
                    tmp = icsTime(200:maxFrame-round(maxFrame/100),m);
                end
                trTmp(:,m) = tmp(301:end);
                trTmp(:,m) = (trTmp(:,m)-mean(trTmp(:,m)))./mean(trTmp(:,m));
                trDisp(:,m) = (trTmp(:,m)-mean(trTmp(:,m)))./std(trTmp(:,m));
                if skewness(trTmp(:,m)) < 0
                    trTmp(:,m) = trTmp(:,m).*-1;
                    trDisp(:,m) = trDisp(:,m).*-1;
                end
            end
            icsTime = trTmp;
            nIcs = size(trace{n},2);
            icsSpace = space{n};
            
            figure(3); clf
            stackplot(trDisp(:,1:min(30, size(icsTime,2))));  % plot up to 30 independent components
            
            figure(4); clf
            for j = 1:size(icsSpace,2)
                subplot(ceil(nIcs/5),5,j);
                imshow2(flipud(toimg(icsSpace(:,j), size(M1,1), size(M1,2))), []);
                title(j)
            end
            keeper = input('keep traces ', 's');
                
            
            if ~isempty(keeper)
                keep = str2num(cell2mat(regexp(keeper,'\d','match')));
                if ~isempty(regexp(keeper,'r'))
                    trace{n}  = icsTime(:,keep).*-1;
                else
                    trace{n}  = icsTime(:,keep);
                end
                space{n} =  icsSpace(:,keep);
            else
                trace{n} = zeros(size(icsTime(:,1)));
                space{n} =  zeros(size(icsSpace(:,1)));
                keep = 0;
            end
            kept{n} = keep;
            
        end
        tr = cat(2,trace{:});
        del = [];
%         tr = ([trace{:}]);
%         ccm = (corrcoef(tr));
%         ccm = triu(ccm,1);
%         [i,j] = ind2sub(size(ccm),find(ccm>0.9));
%         figure(99)
%         p = [i j];
%         for n = 1:length(j)
%             stackplot2(tr(:,p(n,:)))
%             del(n) = input('discard 1 or 2? ');
%             del(n) = p(n,del(n));
%         end
%         if ~isempty(del)
%             tr(:,del) = [];
%         end
        save(strrep(trace_files(k).name, 'trace', 'selected'), 'trace', 'space', 'tr', 'del', 'kept')
    end
end