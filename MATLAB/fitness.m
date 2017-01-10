% fitness(E, P, naa) returns the prediction fitness given energies in
% vector E (naa rows per site, one for each amino acid) and the binary
% observation vector P (a single 1 in each block of naa rows).
function [F, J] = fitness(x, M, nzM, Ec, I, naa, varargin)
persistent Jc weightsPerAA;

usePenalty = 1;
useWeights = 1;

% previosuly use type = 2
np = length(x);             % number of parameters
N = size(M, 1)/naa;         % number of sites (positions)

E = M*x +Ec;
boltzFact = reshape(exp(-E), naa, N);
Q = sum(boltzFact, 1)';
p = reshape(boltzFact ./ repmat(Q', naa, 1), naa*N, 1);

%f = log(boltzFact(I) ./ Q) = log(boltzFact(I)) - log(Q) = -E(I) - log(Q);
f = -E(I) - log(Q);

if (useWeights)
    weights = varargin{1};
    weights = weights/sum(weights);
else
    weights = 1/N;
end
F = -sum(f .* weights);

if (usePenalty)
    % apply regularization; penalize if self terms too far from zero, 
    % and pair terms too far from contact potential
    defaultParams = varargin{2};
    currentParams = varargin{3};
    pUsage = varargin{4};
    optingCols = varargin{5};
    lambda = varargin{6};
%     isself = varargin{7};
    updatedParams = currentParams;
    updatedParams(optingCols) = x;
    penalty = mean((updatedParams - defaultParams) .^ 2 ./ (pUsage + 1)') * lambda;
%     selfpenalty = lambda * sum( (updatedParams(1:20) - defaultParams(1:20)) .^ 2 ./ (pUsage(1:20) + 1)');
%     if length(defaultParams) > 20
%         pairpenalty = 0.025 * lambda * sum( (updatedParams(21:end) - defaultParams(21:end)) .^2 ./ (pUsage(21:end) + 1)');
%     end
%     penalty = selfpenalty + pairpenalty;
    F = F + penalty;
end

% were we asked for the Jacobian?
if (nargout > 1)
    % -- some variable that will not change during optimization
    if (isempty(Jc))
        dEdx = M(I, :); % partial derivative of energies with respect to parameters is just the model matrix
        Jc = zeros(1, np);
        for i = 1:np
            Jc(i) =  sum(dEdx(:, i) .* weights);
        end
        if (useWeights)
            weightsPerAA = weights(reshape(repmat(1:N, naa, 1), naa*N, 1));
        else
            weightsPerAA = 1/N;
        end
    end
    
    % -- finally, partial derivative of frequency of native amino acid at
    % each site
    J = zeros(1, np);
    pw = p .* weightsPerAA;
    for i = 1:np
        J(i) =  Jc(i) - sum(nzM{i}(:, 1) .* pw(nzM{i}(:, 2)));
    end
    
    % Jacobian also changes with the penalty term
    if (usePenalty == 1)
%         if isself
%             J = J + 2 * lambda * (x - defaultParams(optingCols))' ./ (pUsage(optingCols) + 1);
%         else
%             J = J + 0.025 * 2 * lambda * (x - defaultParams(optingCols))' ./ (pUsage(optingCols) + 1);
%         end
        J = J + (x - defaultParams(optingCols))' ./ (pUsage(optingCols) + 1) * 2 * lambda / length(defaultParams);        
    end
end

