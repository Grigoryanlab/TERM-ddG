function [ leftMat, rightVec, paramsUsage ] = loadEnsemble( currentResidue, conNames, header )

self = 20;
pair = 400;

files1 = dir(sprintf('*_%s.pdb', currentResidue));
files2 = dir(sprintf('*_%s_*.pdb', currentResidue));
files = {files1.name, files2.name};

%% lay out the residues coverred by a pdb file
residueLists = {};
for i = 1: length(files)
    fid = fopen(files{i});
    pdbInfo = textscan(fid, '%s', 'Delimiter', '', 'Headerlines', 1);
    fclose(fid);
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

% this part should be in general same
for i = 1: length(files)
    seqf = strcat(header, '_', strrep(files{i}, 'pdb', 'seq'));
    if exist(seqf, 'file') == 0
        continue;
    end
    fileSplit = regexp( strrep(files{i}, '.pdb', ''), '_', 'split');
    cencolid = find(strcmp(residueLists{i}, currentResidue)) + 1;
    
    conInThisFile = [];
    if length(fileSplit) > 2
        importantResidue = fileSplit(2:end);
        conInThisFile = importantResidue(~strcmp(importantResidue, currentResidue));
        conInThisFile = conInThisFile(ismember(conInThisFile, conNames));
        concolid = zeros(size(conInThisFile));
        for c = 1: length(conInThisFile)
            concolid(c) = find(strcmp(residueLists{i}, conInThisFile(c))) + 1;
        end
    end
    
    % test: only use top N sequences;
%     TopN = 20000;
    
    % read the seq files
    nfield = length(residueLists{i}) + 1;
    fid = fopen(seqf);
    seqResults = textscan(fid, repmat('%s ' , [1, nfield]));
    fclose(fid);
    
    nm = length(seqResults{1});
    nc = length(conInThisFile);
    
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
            
            try
                seeds = self + pair * (conid -1) + self * (conaaind -1) + 1;
            catch
                disp(files(i));
            end
            dim2( (d2+1):(nc+1):end ) = reshape(repmat(seeds, 1, self)' + repmat(0:self-1, nm, 1)', nm*self, 1);
        end
    end
    
    smat = sparse(dim1, dim2, ones(size(dim1, 1), 1), nm*self, size(conNames, 1)*pair + self);
    
    % concatenate to existing data
    Sparse = [Sparse; smat];
    rightVec = [rightVec; rmat];
end

leftMat = [];
paramsUsage = [];
%% check the parameters usage, and remove the not used columns
if ~isempty(Sparse)
    paramsUsage = sum(Sparse(logical(rightVec(:,1)), :), 1);
    leftMat = Sparse;
end
%% a side function to deal with non-canonical amino acid
    function [aainds] = aaIndex(names)
        [~, aainds] = ismember(names, AA); 
        aainds(aainds == 21) = 11; % MSE
        aainds(aainds == 0) = 1; % other made ALA
    end

end

