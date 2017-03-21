function [ leftMat, rightVec, defaultParams, paramsUsage, conresidUsed] = loadSearchData( folder, header, conpot )

startedfiles = dir(sprintf('%s/.started.*', folder));
for i = 1:length(startedfiles)
    finishedfile = sprintf('%s/%s', folder, strrep(startedfiles(i).name, 'start', 'finish'));
    if ~exist(finishedfile)
        error('It does not look like every search was finished successfully - %s was found but %s was not, meaning not all jobs in the %s folder reported success', startedfile, finishedfile, folder)
    end
end

maxC = 2;
%cd(folder);
self = 20;
pair = 400;
cennum = regexp(folder, '_', 'split');
cennum = cennum{end};
files = dir(sprintf('%s/*.pdb', folder));

conNames = {};

usefiles = zeros(length(files), 1);
% read the residue id of contacts
icon = 0;
for i = 1:length(files)
    ncon = regexp( strrep(files(i).name, '.pdb', ''), '_', 'split');
    if length(ncon) <= maxC + 2
        usefiles(i) = 1;
    end
    if length(ncon) ~= 3
        continue
    else
        icon = icon + 1;
        conNames{icon} = ncon{3};
    end
end
usedFileInds = find(usefiles == 1);

%% lay out the residues coverred by a pdb file
residueLists = {};
for i = 1: length(usedFileInds)
    pdbInfo = textscan(fopen(sprintf('%s/%s', folder, files(i).name)), '%s', 'Delimiter', '', 'Headerlines', 1);
    pdbInfo = pdbInfo{1};
    residueNums = {};
    for j = 1:length(pdbInfo)
        resinum = strrep(pdbInfo{j}(22:26), ' ', '');
        residueNums = [residueNums, resinum];
    end
    residues = unique(residueNums, 'stable');
    residueLists{i} = residues;
end

AA = {'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR', 'MSE'};

%% iterate over seq files
Sparse = [];
rightVec = [];

for i = 1: length(usedFileInds)
    fInd = usedFileInds(i);
    seqf = strcat(header, '_', strrep(files(fInd).name, 'pdb', 'seq'));
    if exist(sprintf('%s/%s', folder, seqf), 'file') == 0
        continue;
    end
    fileSplit = regexp( strrep(files(fInd).name, '.pdb', ''), '_', 'split');
    cencolid = find(strcmp(residueLists{i}, cennum)) + 1;
    
    if length(fileSplit) > 2
        conInThisFile = fileSplit(3:end);
        concolid = zeros(size(conInThisFile));
        for c = 1: length(conInThisFile)
            concolid(c) = find(strcmp(residueLists{i}, conInThisFile(c))) + 1;
        end
    end
    
    % test: only use top N sequences;
%     TopN = 20000;
    
    % read the seq files
    nfield = length(residueLists{i}) + 1;
    seqResults = textscan(fopen(sprintf('%s/%s', folder, seqf)), repmat('%s ' , [1, nfield]));
%     seqResults = textscan(fopen(sprintf('%s/%s', folder, seqf)), repmat('%s ' , [1, nfield]), TopN);

    nm = length(seqResults{1});
    nc = length(fileSplit)-2;
    
    % create left sparse matrix
    dim1 = reshape( repmat(1:nm*self, nc+1, 1), nm*(nc+1)*self, 1);
    dim2 = zeros(size(dim1));
    
    %% read sequence data
    seqid = 1:nm;
    seqid = seqid';
    cenaa = seqResults{cencolid}(seqid);
    cenaaind = aaIndex(cenaa);
    
    % fills right vector
    rmat = zeros(nm*self, 3);
    rmat(self*(seqid - 1) + cenaaind, 1) = 1;
    rmat(1:nm*self, 2) = ones(nm*self, 1) * i;
    rmat(1:nm*self, 3) = ones(nm*self, 1) * nc;
    
    % fills the self/column part of left matrix
    dim2( 1:(nc+1):end ) = repmat((1:self)', nm, 1);
    
    % contact positions
    if nc > 0
        for d2 = 1:nc
            conaa = seqResults{concolid(d2)}(seqid);
            conaaind = aaIndex(conaa);
            conid = find(strcmp(conNames, conInThisFile(d2)));
            
            seeds = self + pair * (conid -1) + self * (conaaind -1) + 1;
            dim2( (d2+1):(nc+1):end ) = reshape(repmat(seeds, 1, self)' + repmat(0:self-1, nm, 1)', nm*self, 1);
        end
    end
    
    smat = sparse(dim1, dim2, ones(size(dim1, 1), 1), nm*self, icon*pair + self);
    
    % concatenate to existing data
    Sparse = [Sparse; smat];
    rightVec = [rightVec; rmat];
end

% control, only use low level fragments
% usedRows = find(sum(Sparse, 2) <=1);
% % usedRows = 1:size(rightVec, 1);
% Sparse = Sparse(usedRows, :);
% rightVec = rightVec(usedRows, :);

%% check the parameters usage, and remove the not used columns
paramsUsage = sum(Sparse(logical(rightVec(:,1)), :), 1);
leftMat = Sparse;

%% read contact potential approritately
[~, name, ~] = fileparts(folder);
consinfo = textscan(fopen(sprintf('%s/%s.conlist', folder, name)), '%s%s%s%f%s%s');
connames_here = strrep(consinfo{3}, ',', '');
condegree = consinfo{4};
[conInUse, tempind] = ismember(connames_here, conNames);
tempind = tempind(logical(tempind));
condegreeUsed(tempind) = condegree(conInUse);

conres = consinfo{6};
conresid = zeros(length(conres),1);
for c = 1: length(conres)
    conresid(c) = aaIndex(conres{c});
end
conresidUsed(tempind) = conresid(conInUse); 

backgroundConPot = conpotPrior(condegreeUsed, conpot);
defaultParams = [zeros(self, 1); backgroundConPot];

%% a side function to deal with non-canonical amino acid
    function [aainds] = aaIndex(names)
        [~, aainds] = ismember(names, AA); 
        aainds(aainds == 21) = 11; % MSE
        aainds(aainds == 0) = 1; % other made ALA
    end
end

