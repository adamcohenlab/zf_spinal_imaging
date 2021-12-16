function ICA_PCA_array(filenum)
pth = [pwd '/'];
files = dir([pth '*_rois.mat']);
if ~isempty(files)
    %         try
    rois = load([files(filenum).folder '/' files(filenum).name]);
    a = rois.rois;
    if ~exist([files(filenum).folder '/' strrep(files(filenum).name, 'rois', 'ica_pca_trace')], 'file')
        try
            name = strrep(files(filenum).name,'_rois.mat','.nrrd');
            %                     name = strrep(name,'.mat','.cxd');
            disp(['opening image ' pth name]);
            %                     img = bfopen([pth name]);
            %                     dat = img{1,1};
            %                     imgMat = cat(3,dat{:,1});
            %                     totTime = img{1,2}.get(['Global Field ' num2str(size(imgMat,3)) ' Time_From_Start']);
            %                     sr = size(imgMat,3)/totTime;
            imgMat = nrrdread([pth name]);
            %                     maxFrame = floor(30*sr);
            maxFrame = size(imgMat,3);
        catch
            warning(['error loading ' pth name])
        end
        trace = [];
        space = [];
        for n = 1:size(a,2)
            try
            tic
            ymin = round(min(a{1,n}(:,1)));
            ymax = min([round(max(a{1,n}(:,1))) size(imgMat,2)]);
            xmin = round(min(a{1,n}(:,2)));
            xmax = min([round(max(a{1,n}(:,2))) size(imgMat,1)]);
%             M1 = imgMat(xmin:xmax,ymin:ymax,500:maxFrame-round(maxFrame/100));
            M1 = imgMat(xmin:xmax,ymin:ymax,:);
            nframe2 = size(M1, 3);
            
            
            % try some spatial filtering
            
                                M2 = zeros(size(M1));
                                for j = 1:nframe2
                                    M2(:,:,j) = medfilt2(M1(:,:,j), [3 3], 'symmetric');
                                end
            
            
                                movVec = tovec(double(M2));
%             movVec = tovec(double(M1));
            % high-pass filter in time
            
            filtT = 201; % time window for median filtering
            mov2DFilt = movVec - imfilter(movVec, ones(1, filtT)/filtT, 'replicate');
            mov2DFilt = mov2DFilt(:,500:maxFrame-round(maxFrame/100));
            nEigUse = 10; nIcs = 5;
            % PCA
            movVec = double(mov2DFilt);
            covmat = movVec'*movVec;
            % [V, D] = eig(covmat);
            [V, ~] = eigs(covmat,nEigUse);
            V = V(:,end:-1:1);
            u = movVec*V(:,1:nEigUse);
            out2 = double(M1);
            vOrig = tovec(out2 - repmat(mean(out2, 3),[1 1 size(out2,3)]))'*u; 
            % ICA
            try
                [icsTime, ~, sepmat] = sorted_ica(V(:,1:nEigUse),nIcs);
                icsSpace = u(:,1:nEigUse)*sepmat';
                icsTimeOrig = (sepmat*vOrig(:,1:nEigUse)')';
                
            catch
                % if the ICA fails, just give the first pc back, it is
                % probably useless but we don't want an empty matrix as
                % result
                warning('ICA failed');
                icsSpace = u(:,1);
                icsTime = vOrig(:,1);
            end
%             tr = [];
%             for k = 1:size(icsSpace,2)
%                 tr(:,k) = mean(tmpM.*icsSpace(:,k));
%             end
            trace{n}  = icsTimeOrig;
            space{n} =  icsSpace;
            toc
            catch
                warning(['problem on iteration ' num2str(n)])
            end
        end
        disp(['saving ' strrep(files(filenum).name, 'rois', 'ica_pca_trace')])
        save([files(filenum).folder '/' strrep(files(filenum).name, 'rois', 'ica_pca_trace')], 'trace');
        save([files(filenum).folder '/' strrep(files(filenum).name, 'rois', 'ica_pca_space')], 'space');
    end
    %         catch
    %             warning(['problem on iteration ' str2num(filenum)])
    %         end
else
    disp('no files found')
end
exit