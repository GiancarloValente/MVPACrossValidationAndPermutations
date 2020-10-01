dirdata = ' '; % where the downloaded datasets are
dirmatdata = ' '; % where the extracted data will be stored 
filenames = dir(fullfile(dirdata,'*.zip'));

numFiles = numel(filenames);
fileexp  = '_tfMRI_RELATIONAL_level2_hp200_s2_MSMAll.dscalar.nii';


cd('C:\Program Files\7-Zip');% I am using 7-Zip to extract the surface data 
%from the zip files, many other options are available...
for idx = 1:numFiles
    thisfile = fullfile(filenames(idx).folder,filenames(idx).name);
    thissubj = filenames(idx).name(1:6);
    giftiname = [ thissubj  fileexp];
    
    if ~exist(fullfile(dirmatdata,giftiname),'file')
        cmdzip = ['7z.exe e ' thisfile ' -o' dirmatdata ' ' giftiname ' -r > NUL'];
        system(cmdzip);
    end
end
cd(dirmatdata);

%%
% make sure you add the cifti-matlab importer to the path!
% https://github.com/Washington-University/cifti-matlab
% here we used the fieldtrip implementation


EVNames     = {'match','rel'};
ciftinames = dir(fullfile(dirmatdata,'*RELATIONAL*.nii'));
numcifti   = numel(ciftinames);
contrastmapname = @(s1,s2) ['x' s1 '_tfmri_relational_level2_' s2 '_hp200_s2_msmall'];


s0           = ft_read_cifti(fullfile(dirmatdata,ciftinames(1).name));
verticesused = s0.brainstructure<3; % selecting only cortical voxels

numvertices = sum(verticesused);
numcond     = numel(EVNames);
features    = zeros(numel(EVNames)*numcifti,numvertices);
labels      = zeros(numel(EVNames)*numcifti,2);


for idSub    = 1:numcifti
    if mod(idSub,20)==0
        fprintf('.');
    end
    s      = ft_read_cifti(fullfile(dirmatdata,ciftinames(idSub).name));
    thissubj = ciftinames(idSub).name(1:6);
  
    for idEV = 1:numel(EVNames)
            features((idSub-1)*numel(EVNames)+idEV,:) = s.(contrastmapname(thissubj,EVNames{idEV}))(verticesused);
    end
    labels((idSub-1)*numel(EVNames)+[1:numel(EVNames)],1) = [1:numel(EVNames)];
    labels((idSub-1)*numel(EVNames)+[1:numel(EVNames)],2) = idSub;
    
end
fprintf('\n');        

save('dataConnectomeRELATIONAL.mat','features','labels','s0','-v7.3');
% note: this file is large (> 1GB), but it can be recreated from the 889 subjects 
% of the human connectome database (see SubjectIDRELATIONAL.mat)    
  