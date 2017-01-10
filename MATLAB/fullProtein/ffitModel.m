% clear;
function [optParams, obj] = ffitModel(M, rightVec, defaultParams, currentvalue, EAAfreq, paramsUsage, lambda)
naa = 20;
nc = (size(M, 2) - naa)/(naa * naa);
assert(round(nc) == nc, 'could not figure out the number of contacts');
np = size(M, 2);

exclCols = find(paramsUsage == 0); % excluded columns

optParams = currentvalue';

% DerivativeCheck 
% opts = optimoptions('fminunc', 'Algorithm','trust-region', 'Display', 'iter', 'GradObj', 'on', 'DerivativeCheck', 'on');
% opts = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'GradObj', 'on', 'DerivativeCheck', 'on', 'FinDiffRelStep', 10^-6);
opts = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'GradObj', 'on', 'TolFun', 10^-4);

% weigh all queries roughly equally
weights = rightVec(:, 2);
nm = histc(weights, 1:max(weights))/naa;
weights = nm(weights);
weights = 1./(weights + 1000);
notConv = 1;
I = logical(rightVec(:, 1));
weights = reshape(weights, naa, length(weights)/naa);
weights = mean(weights)'; % give a single weight to each position
pobj = inf;

% M is very sparse; it will be useful in objective-function evaluation to
% also separately store its non-zero elements
nzM = cell(1, np);
for i = 1:np
    nz = find(M(:, i));
    nzM{i} = [full(M(nz, i)) nz];
end

%% iteratively optimize
while (notConv)
    disp('optimizing self energies...');
    notConv = 0;
    pFit = 1:naa;
    pFit = setdiff(pFit, exclCols);
    pNotFit = setdiff(1:length(optParams), pFit);
    constE = EAAfreq + M(:, pNotFit) * optParams(pNotFit);
    Mpar = M(:, pFit);
    nzMpar = nzM(pFit);
    if (isempty(pFit)) % extreme case; very sparse data; and if no self parameters, there must be no pair parameters as well
        break;
    end
    % notice defaultParams and paramsUsage; although not all parameters
    % changed on each stage, the penalty of all parameters should be
    % considered on in objective function
    [wSelf, obj] = fminunc(@(x) fitness(x, Mpar, nzMpar, constE, I, naa, ...
        weights, defaultParams, optParams, paramsUsage, pFit, lambda), optParams(pFit), opts);
    clear fitness; % clear persistent variables from inside fitness

    % check if changed
    if (norm(wSelf - optParams(pFit)) > 0.1)
        notConv = 1;
    end
    optParams(pFit) = wSelf;
    
    % optimize contact energies, one contact at a time
    for i = 1:nc
        disp(sprintf('optimizing pair energies for contact %d...', i));
        pFit = naa + 1 + naa*naa*(i-1) : naa + naa*naa*i;
        pFit = setdiff(pFit, exclCols);
        pNotFit = setdiff(1:length(optParams), pFit);
        constE = EAAfreq + M(:, pNotFit) * optParams(pNotFit);
        
        if (isempty(pFit)) % if no pair parameters, move to the next contact
            continue 
        end
        Mpar = M(:, pFit);
        nzMpar = nzM(pFit);

        [wPair, obj] = fminunc(@(x) fitness(x, Mpar, nzMpar, constE, I, naa, ...
            weights, defaultParams, optParams, paramsUsage, pFit, lambda),  optParams(pFit), opts);
        clear fitness; % clear persistent variables from inside fitness

        % check if changed
        if (norm(wPair - optParams(pFit)) > 0.1)
            notConv = 1;
        end
        optParams(pFit) = wPair;
    end
        
    % break if the objective did not improve much
    disp(sprintf('objective from previous cycle is %f, and now is %f', pobj, obj));
    %break % temp
    if (pobj - obj < 0.01)
        break;
    end
    if (~notConv)
        disp('parameters have converged');
    end
    
    pobj = obj;
    
end
end
