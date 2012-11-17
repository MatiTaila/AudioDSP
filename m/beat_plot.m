function beat_plot(b,s,tt,ff,D)
% beat_plot(b,s,tt,ff,D)
%   Plot vertical lines corresponding to the beat times in b on the
%   current plot.  Optional s is the plot style.
%   Optional tt,ff,D define underlying image for imagesc(tt,ff,D)
% 2012-03-27 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; s = ''; end
if nargin < 3; tt = []; end
if nargin < 4; ff = []; end
if nargin < 3; D = []; end

if length(s) == 0; s = '-r'; end

if length(D) > 0
  imagesc(tt,ff,D);
  axis('xy');
  colormap(1-gray);
end

ax = axis;

hold on;
plot([b;b], repmat([ax(3) ax(4)]',1,length(b)), s);
hold off;
