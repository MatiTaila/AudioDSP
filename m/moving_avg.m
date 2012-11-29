function y = moving_avg(x, avg_size)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% function y = moving_avg(x, avg_size)
%
% Applies a moving average filter to the input signal.
%
% Inputs:
%   x: Input signal.
%   avg_size: Number of samples to use in the moving average.
% Outputs:
%   y: Filtered signal.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

y = zeros(length(x),1);
for i=avg_size+1:length(x)-avg_size
    frame = x(i-avg_size:i+avg_size);
    y(i) = mean(frame);
end
y(1:avg_size) = ones(avg_size,1)*y(avg_size+1);
y(end-avg_size+1:end) = ones(avg_size,1)*y(end-avg_size);