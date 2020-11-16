function [x,l]                  = generatedatamultirun(Params)
% code to generate simulated data according to the values provided in Params
%
% Copyright (c) Giancarlo Valente 2020
% giancarlo.valente@maastrichtuniversity.nl
%
% Giancarlo Valente licenses this file to you under the MIT License.
% See the LICENSE file for more information

nruns                           = Params.nRuns;
nsamplesrunClass1               = floor(Params.nSamplesPerClass1./nruns);
nsamplesrunClass2               = floor(Params.nSamplesPerClass2./nruns);

if (nsamplesrunClass1*nruns     ~= Params.nSamplesPerClass1) || ...
        (nsamplesrunClass2*nruns   ~= Params.nSamplesPerClass2)
    disp('Warning: Given the number of runs and the number of samples per class, not all the runs will have the same examples');
    disp('Warning: The last run will have more samples than the others');
end

l1                              = zeros(Params.nSamplesPerClass1,1); % l1 contains the labels of run membership for class1
xClass1                         = zeros(Params.nSamplesPerClass1,Params.nVoxels); % xClass1 contains the data of class1
l2                              = zeros(Params.nSamplesPerClass2,1); % l1 contains the labels of run membership for class2
xClass2                         = zeros(Params.nSamplesPerClass2,Params.nVoxels); % xClass1 contains the data of class2

for irun                        = 1:nruns
    if irun                     == nruns
        k1                       = nsamplesrunClass1*(irun-1)+1   : Params.nSamplesPerClass1;
        k2                       = nsamplesrunClass2*(irun-1)+1   : Params.nSamplesPerClass2;
                
    else
        k1                       = nsamplesrunClass1*(irun-1)+1 : nsamplesrunClass1*irun;
        k2                       = nsamplesrunClass2*(irun-1)+1 : nsamplesrunClass2*irun;
                
    end
    thisRunVariance              =  exp(randn*Params.runVariance); % using a log-normal distribution to model run-specific variance
    if strcmpi(Params.distributionClass1.DistributionName,'normal')&& Params.runVariance>0
        Params.distributionClass1.Parameters{2} =thisRunVariance;
    end
    if strcmpi(Params.distributionClass2.DistributionName,'normal')&& Params.runVariance>0
        Params.distributionClass2.Parameters{2} = thisRunVariance;
    end
    
    % generating data according to the distribution provided in Params
    cellargument                = [{Params.distributionClass1.DistributionName} ...
        Params.distributionClass1.Parameters {[numel(k1) Params.nVoxels]}];
    xClass1(k1,:)        = random(cellargument{:});
    
    % 'injecting a unviariate difference in the discriminative voxels
    % (works only in new Matlab versions, otherwise use repmat or bsxfun!
    xClass1(k1,:)        = xClass1(k1,:) +  [ones(1,Params.nDiscriminativeVoxels)*Params.differenceMagnitude/2 ,...
        zeros(1,Params.nVoxels-Params.nDiscriminativeVoxels)];
    l1(k1)               = irun;
    
    
    
        % generating data according to the distribution provided in Params
    cellargument                = [{Params.distributionClass2.DistributionName} ...
        Params.distributionClass2.Parameters {[numel(k2) Params.nVoxels]}];
    xClass2(k2,:)        = random(cellargument{:});
    
    % 'injecting a unviariate difference in the discriminative voxels
    % (works only in new Matlab versions, otherwise use repmat or bsxfun!
    xClass2(k1,:)        = xClass2(k1,:) +  [-ones(1,Params.nDiscriminativeVoxels)*Params.differenceMagnitude/2 ,...
        zeros(1,Params.nVoxels-Params.nDiscriminativeVoxels)];
    l2(k1)               = irun;
    
    
    
end

x                               = [xClass1;xClass2];
l                               = [ [ones(size(l1)) l1]; [2*ones(size(l2)) l2]];

