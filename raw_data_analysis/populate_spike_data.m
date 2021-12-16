 for n = 1:length(data)
    if ~isempty(data(n).av)
        tmp = findSpikes(data(n).trace);
        data(n).spikes = tmp;
    end
end
save('complete_dataset.mat', 'data')