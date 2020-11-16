function [l] = generateLabels(Params)
% generate labels with three designs: random, blocked or alternated
%
% Copyright (c) Giancarlo Valente 2020
% giancarlo.valente@maastrichtuniversity.nl
%
% Giancarlo Valente licenses this file to you under the MIT License.
% See the LICENSE file for more information

l   = ones(Params.nSamplesperRun,Params.nRuns);

switch lower(Params.designType)
    case 'random'
        
        l(Params.nSamplesperRun/2+1 : end,:) = 2;
        for idx = 1:Params.nRuns
            
            l(:,idx) = l(randperm(Params.nSamplesperRun),idx);
        end
    case 'alternate'
        l(2:2:end,:) = 2;
        tmp          =binornd(1,.5,[1 size(l,2)]);
        l(:,logical(tmp)) = 3-l(:,logical(tmp));
        
    case 'blocked'
        nblocks  = 2;
        blocksize = round(Params.nSamplesperRun/nblocks);
        tmp = [1:blocksize]+ [0:2:(nblocks-1)]'*blocksize;
        tmp = tmp(:);
        l(tmp,:)=2;
        
        tmp          =binornd(1,.5,[1 size(l,2)]);
        l(:,logical(tmp)) = 3-l(:,logical(tmp));
end

l = [ l(:) kron([1:Params.nRuns]',ones(Params.nSamplesperRun,1))];


    


