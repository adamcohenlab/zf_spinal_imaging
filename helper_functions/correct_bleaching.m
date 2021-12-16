%% correct bleaching in fluorescent signal
% corrects bleaching by creating an envelop around the minimal values of a
% fluorescent signal and fitting a double exponential
% input:
% F_raw = vector with fluorescent signal
% output:
% F_corr = vector of corrected signal

function F_corr = correct_bleaching(F_raw)

% smooth signal to reduce outliers
f1sm=smooth(F_raw,10);

%calculate envelope 
ind=1;
% start at timepoint 1
cc = [];
while ind<length(f1sm)                          % iterate through the signal
    for i=ind:(length(f1sm)-1)                  % calculate all the slopes between the current point and all the future points
        tg(i)=(f1sm(i+1)-f1sm(ind))/(i+1-ind);
    end
    lind=ind;
    [mn,ind] = min(tg(ind:end));                % get the point with the minimal slope
    if mn>0                                     % if the slope is poisitive, use the last one and end iteration
        cc(lind:length(f1sm))=f1sm(lind);
        break
    end
    %clear tg
    ind=ind+lind;                               % the index of the newly found point
    % cc(ind)=f1(ind);
    y1=f1sm(ind);
    y2=f1sm(lind);
    cc(lind:ind)=y2:(y1-y2)/(length(f1sm(lind:ind))-1):y1;  % create a straight line between the last and the new point
end

%fit double exponetial to the envelope
x = (1:length(cc))';
fobj = fit(x,cc','exp2');
%extract final correction curve
cc2 = feval(fobj,x)';


% %check whether fit was ok
% while gof.rsquare<0.9
%         fig = figure;
%         plot(F_raw)
%         hold all
%         plot(fobj)
%         w = zeros(length(F_raw),1);
%         wrange = input('the rsquare value was (below 0.9) do you want to exclude a region from the fit?\n (enter range to exclude or hit return) ');
%         if isempty(wrange)
%             close(fig)
%             break
%         elseif wrange>length(F_raw)
%             wrange = input('out of range, enter new range ');
%         end
%         w(wrange) = 1;
%         [fobj, gof] = fit(x,F_raw,'exp2','Exclude',w);
%         close(fig)
%         %better = input('better? (y/n) ');
% end

%extract curve from fitting data




% correct signal by adding difference between decay and starting value
F_corr = F_raw + (cc2(1)-cc2)';


%F_corr3 = F_raw./corr_exitcurve2;
% figure
% plot(F_raw)
% hold all
% plot(f1sm)
% plot(cc)
% plot(cc2)
% plot(F_corr)
end