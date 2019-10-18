function Ylabels = mapToLabels( Y, uniqueLabels )

Ylabels = zeros(length(Y),length(uniqueLabels));

for k = 1:length(uniqueLabels)
    Ylabels(:,k) = abs(Y-uniqueLabels(k));
end

[~,Ylabels] = min(Ylabels,[],2);

Ylabels = uniqueLabels(Ylabels);

return