function [X_train,Y_train,X_test,Y_test,Names] = pumadyn_CreateDataSets

MAINDIR = {'~/Git/DataSets/','/Volumes/Large Storage/DataSets/pumadyn/','/Volumes/FastStorage/Data Sets/pumadyn/','/ztmp/DataSets/pumadyn'};

X_train   = cell(8,1);
Y_train   = cell(8,1);
X_test    = cell(8,1);
Y_test    = cell(8,1);
Names     = cell(8,1);

for k = 1:length(MAINDIR)
    lDirs = dir([MAINDIR{k} 'pumadyn-*']);
    if ~isempty(lDirs), break; end
end
if isempty(lDirs), warning('\n pumadyn_CreateDataSets: Could not find data directory. Please add it to this function.'); return; end;
mainDir = MAINDIR{k};

for k = 1:length(lDirs)
    Z = load([mainDir lDirs(k).name '/Dataset.data']);
    X_train{k} = Z(:,1:end-1)';
    Y_train{k} = Z(:,end);
    Z = load([mainDir lDirs(k).name '/accel/Prototask.data']);
    X_test{k} = Z(:,1:end-1)';
    Y_test{k} = Z(:,end);    
    Names{k} = lDirs(k).name;
end

return