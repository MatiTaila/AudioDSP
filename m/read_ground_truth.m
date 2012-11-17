function [T,tempos] = read_ground_truth(file)
 D = textread(file);
 tempos = zeros(size(D,1),1);
for i =1:size(D,1)
    T{i} = D(i,D(i,:)~=0);
    tempos(i) = median(diff(T{i}));
end
