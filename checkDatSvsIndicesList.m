function [isMatch, indicesList, fileList] = checkDatSvsIndicesList(matrix, tolerance, datFile)
% CHECKDATSVSINDICESLIST  Compare groups stored in <datFile> with
%                        the groups returned by findZerosAndSimilarEntries.
%
%   isMatch = CHECKDATSVSINDICESLIST(M, tol, 'chsh.dat-s')
%
% * Works with arbitrary headers, blank lines and extra whitespace.
% * Reading is based on a regular expression – no silent truncation.
% * The comparison is insensitive to the ordering of groups or of
%   indices inside each group.

% -------------------------------------------------------------------------
% ❶  Obtain indicesList produced by YOUR routine
% -------------------------------------------------------------------------
[indicesList, ~] = findZerosAndSimilarEntries(matrix, tolerance);

% -------------------------------------------------------------------------
% ❷  Parse the .dat-s file  →  numeric triplets  [ ID  row  col ]
% -------------------------------------------------------------------------
fid = fopen(datFile, 'rt');
if fid < 0
    error('Cannot open "%s".', datFile);
end
rawTxt = fread(fid, '*char')';
fclose(fid);

% ― one regexp to grab every line that *begins* with an integer ------------
expr      = '(?m)^\s*(\d+)\s+\d+\s+(\d+)\s+(\d+)';   % ID blk row col …
tokens    = regexp(rawTxt, expr, 'tokens');

if isempty(tokens)
    error('No numeric entries recognised in "%s".', datFile);
end

trip      = str2double( vertcat(tokens{:}) );        % [ID  row  col]
IDs       = trip(:,1);
rowsFile  = trip(:,2);
colsFile  = trip(:,3);

% -------------------------------------------------------------------------
% ❸  Build fileList: one cell per unique ID (867 for the CHSH test file)
% -------------------------------------------------------------------------
uIDs      = unique(IDs,'stable');
fileList  = cell(numel(uIDs),1);
for k = 1:numel(uIDs)
    msk          = IDs == uIDs(k);
    fileList{k}  = [ rowsFile(msk) , colsFile(msk) ];
end

% -------------------------------------------------------------------------
% ❹  Canonicalise & compare ------------------------------------------------
% -------------------------------------------------------------------------
indicesList = canonicalise(indicesList(:));   % force column vector of cells
fileList    = canonicalise(fileList(:));

isMatch = isequal(indicesList, fileList);

if isMatch
    fprintf('✅  Lists match (%d groups).\n', numel(indicesList));
else
    reportMismatch(indicesList, fileList);
end
end
% =========================================================================
function L = canonicalise(Lin)
% Sort every [row col] matrix, then lexicographically sort the cells.
for k = 1:numel(Lin)
    Lin{k} = sortrows(Lin{k}, [1 2]);
end
keys     = cellfun(@(m) sprintf('%d,%d;', m.'), Lin, 'uni', false);
[~, ord] = sort(keys);
L        = Lin(ord);
end
% -------------------------------------------------------------------------
function reportMismatch(A, B)
% Print a compact diagnostic summary.
extraA   = setdiff(cellfun(@dump, A, 'uni', false), ...
                   cellfun(@dump, B, 'uni', false));
extraB   = setdiff(cellfun(@dump, B, 'uni', false), ...
                   cellfun(@dump, A, 'uni', false));

fprintf('❌  Lists differ.\n');
fprintf('   Groups in indicesList but not in *.dat-s*: %d\n', numel(extraA));
fprintf('   Groups in *.dat-s*  but not in indicesList: %d\n', numel(extraB));

if ~isempty(extraA)
    fprintf('   First missing group (indicesList → dat-s):\n');
    disp( str2num(extraA{1}) ); %#ok<ST2NM>
end
if ~isempty(extraB)
    fprintf('   First extra group (dat-s → indicesList):\n');
    disp( str2num(extraB{1}) ); %#ok<ST2NM>
end
end
% -------------------------------------------------------------------------
function s = dump(m)
% Helper: unique textual signature for a matrix.
s = sprintf('%d,%d;', m.');
end
