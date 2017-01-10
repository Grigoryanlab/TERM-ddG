function [optParams] = main(dirname, header, outFile, lambda, noprior, parUse, contactPotential)

naa = 20;

if isempty(lambda)
    lambda = 1000;
end

tstart = cputime;

%% read input
disp(strcat('start ', dirname));

conpot = '/home/anthill/fzheng/home/Thesis/Data/searchDB/contactPotential/conpot.free.dun.list';
if ~isempty(contactPotential)
    conpot = contactPotential;
end

[leftMat, rightVec, defaultParams, paramsUsage, ~] = loadSearchData(dirname, header, conpot);

%% if don't use contact potential prior
if ~isempty(noprior)
    defaultParams = zeros(size(defaultParams));
end

%% if only cares the parameter usage
if ~isempty(parUse)
    paramsUsage = full(paramsUsage);
    save(strcat(outFile, '.paruse'), 'paramsUsage');
    return;
end

%% prior background values;
EAAfreq = importdata('aafreq.txt');
EAAfreq = -log(EAAfreq);
EAAfreq = EAAfreq - mean(EAAfreq(:));
EAAfreqVec = repmat(EAAfreq, size(rightVec, 1)/naa, 1);

%% optimization
optParams = fitModel(leftMat, rightVec, defaultParams, EAAfreqVec, outFile, paramsUsage, lambda);

disp(strcat('finish ', dirname));
save(strcat(outFile, '.mat'), 'optParams');

telapse = cputime - tstart;
disp(telapse);
