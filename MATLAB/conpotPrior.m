function [ bconpot ] = conpotPrior(condegree, conpot)

naa =20;
nc = length(condegree);

conbin = [0.02 0.05 0.1 0.2 0.5 1];
conpotfiles = textscan(fopen(conpot), '%s');
conpots = cell(5, 1);
for i = 1:5
    conpots{i} = importdata(conpotfiles{1}{i});
end

bconpot = zeros(naa * naa * nc, 1);
bins = zeros(nc, 1);
for i  = 1:nc
    bins(i) = find(conbin < condegree(i), 1, 'last' );
    bconpot(1 + naa*naa*(i-1) : naa*naa*i) = reshape(conpots{bins(i)}, 400, 1);
end

end

