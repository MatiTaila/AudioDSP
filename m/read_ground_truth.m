function [T,tempos,beats] = read_ground_truth(file)
 D = textread(file);
 tempos = zeros(size(D,1),1);
for i =1:size(D,1)
    T{i} = D(i,D(i,:)~=0);
    tempos(i) = median(diff(T{i}));
end

beats = zeros(size(D,2),1);
for i =1:size(D,2)
    beats(i) = median(D(D(:,i)~=0,i));
end

fileID = fopen (['./beat_examples/train/Dataset1/' file(1:end-4) '.lab'],'w');
fprintf (fileID,'%6.2f',beats');
fclose (fileID);