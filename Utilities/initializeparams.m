function Params = initializeparams(varargin)


% setting defaults

Params.nDatasets                        = 30;  % nDatasets are generated
Params.nVoxels                          = 100; % number of voxels (total)
Params.nSamplesPerClass1                = 60; % trials per class 1 (total)
Params.nSamplesPerClass2                = 60; % trials per class 2 (total)
Params.distributionClass1               = struct('DistributionName','normal','Parameters',{{0 1}}); % distribution of voxel values in class 1
Params.distributionClass2               = struct('DistributionName','normal','Parameters',{{0 1}}); % distribution of voxel values in class 2
Params.nRuns                            = 10; % number of acquisition runs. We assume that trials are balanced and equally distributed across runs                                 
Params.nDiscriminativeVoxels            = 20; % voxels exhibiting a univariate difference
Params.differenceMagnitude              = .5; % univariate difference between classes in the discrimanative voxels
Params.runVariance                      = 0; % variance across runs can vary (exp(rand*Paras.runVariance));
Params.crossValidation                  = 'kfold'; % could be kfold, hold out, leave-run-out
Params.nFolds                           = 5; % used if cross-validation is enabled 
Params.percentageHoldOut                = .2; % used if hold out is enabled
Params.repeatCrossValidation            = true; % if meaningful, the training-test division in CV is repeated multiple times
Params.nRepetitionsCrossValidation      = 5; % used if repeatCrossValidation is enabled
Params.repeatHoldout                    = false; % the training test division in hold out is repeated multiple times
Params.nRepetitionsHoldout              = 30; % used if repeatHoldout is true;
Params.algorithm                        = 'svm'; % learning algorithm, could be svm' or 'liblinear' or 'gnb' (parallelized)
Params.nPerm                            = 1e3; % total number of permutations
Params.permutationType                  = 'beforesplit'; % can be 'beforesplit','withinsplit','witinsplittrain','withinsplittest': 
%                     1) beforesplit: first permute the labels, and then test in cross-validation
%                     2) withinsplit: within each split, randomize training and testing separately
%                     3) withinsplittrain: within each split, randomize only training data
%                     4) withinsplittest: within each split, randomize only test data


% checking for modifications of the default parameters
assert(mod(numel(varargin),2) == 0,'Please provide parameters in pairs (Name - Value)');
for ind = 1:2:numel(varargin)
    assert(isfield(Params,varargin{ind}),'Please provide a valid parameter name');
    Params.(varargin{ind}) = varargin{ind+1};
end
    
