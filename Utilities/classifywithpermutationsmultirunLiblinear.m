function [err,eP,varargout] = classifywithpermutationsmultirunLiblinear(x,l,Splits,Params)
% Code to classify using LibLinear 
%
% Copyright (c) Giancarlo Valente 2020
% giancarlo.valente@maastrichtuniversity.nl
%
% Giancarlo Valente licenses this file to you under the MIT License.
% See the LICENSE file for more information


eP                          = [];
eP1                         = [];
eP2                         = [];

xsparse                     = sparse(x);

% running classification
[err,err1,err2]             = classifyLibLinear(xsparse,Splits,l);

% running permutations if needed
if Params.nPerm > 0
    
    eP                      = zeros(Params.nPerm,numel(Splits));
    eP1                     = zeros(Params.nPerm,numel(Splits));
    eP2                     = zeros(Params.nPerm,numel(Splits));
    parfor iPerm            = 1:Params.nPerm
    
        % permuting the labels according to the desired scheme
        lp                  = permuteLabels(l,Params,Splits);
        % evaluating classifier performance
        [terr,terr1,terr2]  = classifyLibLinear(xsparse,Splits,lp);
    
        eP(iPerm,:)         = terr;
        eP1(iPerm,:)        = terr1;
        eP2(iPerm,:)        = terr2;
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









function [e,e1,e2]          = classifyLibLinear(x,Splits,l)


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
    
    mdl                     = train(ltrain, xtrain ,'-s 0 -q');
    pred                    = predict(ltest,xtest,mdl,'-q');
    
    e(iSplit)               = sum(pred ~= ltest);
    e1(iSplit)              = sum(pred(ltest==1) ~= ltest(ltest==1));
    e2(iSplit)              = sum(pred(ltest==2) ~= ltest(ltest==2));
    
end