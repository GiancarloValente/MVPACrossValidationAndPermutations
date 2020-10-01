function [err,eP,varargout] = classifywithpermutationsmultirunLIBSVM(x,l,Splits,Params)

eP                          = [];
eP1                         = [];
eP2                         = [];


% precomputing kernel 
kernel                      = x*x';


% running classification
[err,err1,err2]             = classifyLIBSVM(kernel,Splits,l);

% running permutations if needed
if Params.nPerm > 0
    
    eP                      = zeros(Params.nPerm,numel(Splits));
    eP1                     = zeros(Params.nPerm,numel(Splits));
    eP2                     = zeros(Params.nPerm,numel(Splits));
    parfor iPerm               = 1:Params.nPerm
    
        % permuting the labels according to the desired scheme
        lp                  = permuteLabels(l,Params,Splits);
        % evaluating classifier performance
        [terr,terr1,terr2]  = classifyLIBSVM(kernel,Splits,lp);
    
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









function [e,e1,e2,varargout]= classifyLIBSVM(kernel,Splits,l)


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
    kerneltrain             = [(1:numel(ltrain))' kernel(Splits(iSplit).train,Splits(iSplit).train)];
    kerneltest              = [(1:numel(ltest))' kernel(Splits(iSplit).test,Splits(iSplit).train)];
    
    mdl                     = svmtrain(ltrain, kerneltrain, '-q -t 4'); %#ok<*SVMTRAIN>
    pred                    = svmpredict(ltest, kerneltest, mdl,'-q');
    
    e(iSplit)               = sum(pred ~= ltest);
    e1(iSplit)              = sum(pred(ltest==1) ~= ltest(ltest==1));
    e2(iSplit)              = sum(pred(ltest==2) ~= ltest(ltest==2));
    
end