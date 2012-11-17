function [T,tempo] = beat_ground_truth(F, VERBOSE)
% [T,tempo] = get_ground_truth(F, VERBOSE)
%    F is the name of a mirex06-tyle ground truth file.
%    Return T as a cell array of ground-truth tap times
%    .. but only return tracks whose tempo is within 10 BPM of
%    most popular subjective tempo.
% 2012-03-27 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; VERBOSE = 0; end

if nargout == 0; VERBOSE = 1; end

% File is rectangular array with each row the tap times, padded
% with zeros
D = textread(F);
nD = size(D,1);
tempos = zeros(nD,1);

for i = 1:nD
  db = D(i,:);
  % remove zeros
  db = db(db>0);
  % figure median IOI
  tempos(i) = 60/median(diff(db));
  T{i} = db;
end

% figure consensus tempo
binwidth = 10;
tempi = 10:binwidth:240;
tempohist = hist(tempos,tempi);
% sum up adjacent bins
tempohist = [0,tempohist] + [tempohist,0];
% find largest
[vv,xx] = max(tempohist);
tpmin = tempi(xx-1)-binwidth/2;
tpmax = tempi(xx)+binwidth/2;
% average actual values in those two bins
goodtempos = find( (tempos >= tpmin) & (tempos <= tpmax) );
tempo = mean(tempos(goodtempos));

if VERBOSE
  disp(['Consensus tempo = ',num2str(tempo), ...
        ' BPM (', num2str(length(goodtempos)),'/', ...
        num2str(length(tempos)),' examples)']);
end

% Only return the consistent tempos
T = T(goodtempos);
