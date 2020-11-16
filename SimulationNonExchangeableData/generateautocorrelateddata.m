function [X,runidx]  = generateautocorrelateddata(Params,varargin)
% generate data where consecutive examples have a given correlationS
%
% Copyright (c) Giancarlo Valente 2020
% giancarlo.valente@maastrichtuniversity.nl
%
% Giancarlo Valente licenses this file to you under the MIT License.
% See the LICENSE file for more information


% creating the covariance matrix, with bands
corr_shape = normpdf(-Params.correlationWidth:Params.correlationWidth,0,Params.correlationS);
corr_shape = corr_shape./max(corr_shape);
C = full(spdiags([ones(Params.nSamplesperRun,1)*corr_shape],-Params.correlationWidth:Params.correlationWidth,...
    Params.nSamplesperRun,Params.nSamplesperRun));

if false
    subplot(3,1,3), plot(-Params.correlationWidth:Params.correlationWidth,corr_shape);
    subplot(3,1,1:2); imagesc(C); axis square
end

Y = mvnrnd(zeros(Params.nSamplesperRun,1),C,Params.nVoxels*Params.nRuns)';
tmpY = kron([1:Params.nRuns],ones(Params.nSamplesperRun,Params.nVoxels));
X = zeros(Params.nSamplesperRun*Params.nRuns,Params.nVoxels);
tmpX = kron([1:Params.nRuns]',ones(Params.nSamplesperRun,1))*ones(1,Params.nVoxels);
for idx = 1:Params.nRuns
    X(tmpX==idx) = Y(tmpY==idx);
end
runidx = tmpX(:,1);