clear variables
taskname            = 'relational'; % could be either motor or relational
addmexlibraries;
disp('loading data');
namedataset = ['dataConnectome', upper(taskname) '.mat'];

load(fullfile(pwd,namedataset),'features','labels','s0');
switch lower(taskname)
    case 'motor'
        EVNames         = {'lf','lh','rh','rf','t'};
    case 'relational'
        EVNames         = {'match','rel'};
end
numConds        = numel(EVNames);
switch lower(taskname) 
    case 'motor'
        class1                          = 'rh'; % right hand
        class2                          = 'rf'; % right foot
    case 'relational'
        class1                          = 'match';
        class2                          = 'rel';
end
numiterations                   = 100;  % the whole analysis is repeated 100 times with random selection of subjects and voxels
numpermutations                 = 1000; % for each decoding, 1000 permutations are used]
numvoxselected                  = [10 100 1000];
disp('done');
%%

res(numel(numvoxselected),numiterations) = struct('E',[],'Er',[],'p',[],'prep',[]);
for idvox           = 1:numel(numvoxselected)
    for id              = 1:numiterations
        [res(idvox,id).E,res(idvox,id).Er,res(idvox,id).p,res(idvox,id).prep] ...
            = comparexvalconnectome_iteration(features,labels,EVNames,class1,class2,numpermutations,numvoxselected(idvox));
    end
    namesave = sprintf('Results%s_%dperm_%diter.mat',upper(taskname),numpermutations,numiterations);
    save(namesave,'res','numvoxselected');
end

