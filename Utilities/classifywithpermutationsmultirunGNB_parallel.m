function [err,eP,varargout] = classifywithpermutationsmultirunGNB_parallel(x,l,Splits,Params)
% parallelized implementation of GNB permutations based on Ontivero-Ortega et al. Neuroimage 2017. Please note that the permutations are correct only if the design is balanced within each run.
%
% Copyright (c) Giancarlo Valente 2020
% giancarlo.valente@maastrichtuniversity.nl
%
% Giancarlo Valente licenses this file to you under the MIT License.
% See the LICENSE file for more information


eP                          = [];
eP1                         = [];
eP2                         = [];

% running classification
[err,err1,err2]             = classifyGNB(x,Splits,l);

% running permutations if needed
if Params.nPerm > 0
    
    eP                      = zeros(Params.nPerm,numel(Splits));
    eP1                     = zeros(Params.nPerm,numel(Splits));
    eP2                     = zeros(Params.nPerm,numel(Splits));
    l1                      = l;
    l1(l(:,1) ==2,1)        = 0;
    
   
    for iPerm               = Params.nPerm:-1:1
        lp(:,:,:,iPerm)     = permuteLabels(l1,Params,Splits);
    end
    
    for iSplit              = 1:numel(Splits)
        ltrain              = squeeze(lp(Splits(iSplit).train,1,iSplit,:));
        ltest               = squeeze(lp(Splits(iSplit).test,1,iSplit,:));
        xtrain              = x(Splits(iSplit).train,:);
        xtest               = x(Splits(iSplit).test,:);
        pred                = permutations4GNB(xtest,xtrain,ltrain);
        
        eP(:,iSplit)        = sum(pred ~= ltest)';
        eP1(:,iSplit)       = sum(pred(ltest==1) ~= ltest(ltest==1))';
        eP2(:,iSplit)       = sum(pred(ltest==2) ~= ltest(ltest==2))';
    
    end
end   
    
    

if nargout-2>=1
    varargout{1}            = err1;
end
if nargout-2>=2
    varargout{2}            = err2;
end
if nargout-2 >= 3
    varargout{3}            = eP1;
end
if nargout-2 >= 4
    varargout{4}            = eP2;
end




function [e,e1,e2]          = classifyGNB(x,Splits,l)


e                           = zeros(1,numel(Splits));
e1                          = zeros(1,numel(Splits));
e2                          = zeros(1,numel(Splits));

isPermutation               = false;
if size(l,3)                > 1 
    isPermutation           = true;
end

for iSplit                  = 1:numel(Splits)
    
    if isPermutation
        ltrain              = l(Splits(iSplit).train,1,iSplit);
        ltest               = l(Splits(iSplit).test,1,iSplit);
    else
        ltrain              = l(Splits(iSplit).train,1);
        ltest               = l(Splits(iSplit).test,1);
    end
    xtrain                  = x(Splits(iSplit).train,:);
    xtest                   = x(Splits(iSplit).test,:); 
    
    pred                    = predictgnb_sharedvariance(xtrain,ltrain,xtest);
    
    e(iSplit)               = sum(pred ~= ltest);
    e1(iSplit)              = sum(pred(ltest==1) ~= ltest(ltest==1));
    e2(iSplit)              = sum(pred(ltest==2) ~= ltest(ltest==2));
    
end


function [pred] = predictgnb_sharedvariance(xtrain,ltrain,xtest)

classeslabels = unique(ltrain);
mu1           = mean(xtrain(ltrain==classeslabels(1),:));
mu2           = mean(xtrain(ltrain==classeslabels(2),:));
stdshared     = std([[xtrain(ltrain==classeslabels(1),:) - mu1]; [xtrain(ltrain==classeslabels(2),:) - mu2]]);

logLik1       = sum(-(xtest-mu1).^2./(2.*stdshared),2);
logLik2       = sum(-(xtest-mu2).^2./(2.*stdshared),2);

pred          = nan(size(xtest,1),1);
pred(logLik1 > logLik2) = classeslabels(1);
pred(logLik1 < logLik2) = classeslabels(2);

function C = permutations4GNB(Xtest,X,Y)
% Xtest ntest x p test data
% X     n     x p training data
% Y     n     x nperm Y labels for each permutation
% C     ntest x nperm classification for the ntest samples


%%%%%%%%%%%%%  NOTE: this is only valid only for balanced designs!!!!!!!

n     = size(X,1);
mu    = 2.*mean(X);             % global mean
muX2  = mean(X.^2);             % global mean for X2
X     = 2*X./n;                 % valid only for balanced designs 
muA   = Y'*X;                   % mean  class A
muB   = mu - muA;               % mean  class B
vX    = muX2 - (muA.^2  + muB.^2)/2; %valid only for balanced designs 
% discriminant function linear term
l1   = (muA-muB)./vX;
K   = sum( (    muA.^2 - muB.^2)./vX ,2)./2; % K    = sum( (    muA.^2 - muB.^2)./(2.*vX) ,2); % discriminant constant term
C    = 1.*(Xtest*l1' > K'); % classification

