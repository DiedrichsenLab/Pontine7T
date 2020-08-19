function [R2,R]=bsp_patternConsistency_crossval(Y,partition,condition,varargin)
% function R_cross=rsa_patternConsistency(Y,partition,condition,varargin);
% Calcylates measures of pattern consistency: R2 (SS_between/SS_total) and R (correlation)
% INPUT:
%  Y                : Noise-normalized activation patterns, a KR x P matrix
%  partition        : KR x 1 integer value that indicates the partition for crossvalidation (typically run number), zeros will be ignored
%  condition        : KR x 1 vector of conditions, zeros will be ignored
%
% OUTPUT:
%   R2              : crossvalidated R2 (1-SSR/SST)
%   R               : crossvalidated R (correlation)

% EBerlot, May 2018

removeMean=1; % by default removeMean
vararginoptions(varargin,{'removeMean'});

if ~isnan(Y) % if no data provided
    
    part    = unique(partition)';
    part    = part(part~=0)';
    numPart = numel(part);
    cond    = unique(condition)';
    cond    = cond(cond~=0)';
    numCond = numel(cond);
    
    [~,numVox]   = size(Y);
    
    A = zeros(numCond,numVox,numPart);
    X = indicatorMatrix('identity_p',condition);
    
    % Estimate condition means within each run
    for i=1:numPart
        Xa = X(partition==part(i),:);
        Ya = Y(partition==part(i),:);
        A(:,:,i) = pinv(Xa)*Ya;
        if (removeMean) % per partition
            A(:,:,i)=bsxfun(@minus,A(:,:,i),sum(A(:,:,i),1)/numCond);
        end
    end
    Pred = zeros(size(A));
    % prep betas - leave one out manner
    for i=1:numPart
        trainRun     = part~=i;
        Pred(:,:,i)  = mean(A(:,:,trainRun),3);
    end 
    SSR      = sum(sum(sum((A-Pred).^2)));
    SST      = sum(sum(sum(A.^2)));
    SSP      = sum(sum(sum(Pred.^2))); 
    SSPA     = sum(sum(sum(Pred.*A))); 
    
    R2 = 1-SSR/SST;
    R  = SSPA / sqrt(SST*SSP);
else
    warning('no data given');
    R2 = NaN;
    R = NaN;
end


