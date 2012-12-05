close all
clear all
clc

agents(1).Pm = 100; agents(1).Phi = 134; agents(1).age = 1;
agents(2).Pm = 110; agents(2).Phi = 134; agents(2).age = 2;
agents(3).Pm = 105; agents(3).Phi = 135; agents(3).age = 3;

fs = 44100;
n_hop = 100;
KILLED_BY_REDUNDANCY = 0;
REDUNDANCY_P_MAX     = 11.6e-3*(fs/n_hop);
REDUNDANCY_PHI_MAX   = 23.2e-3*(fs/n_hop);

aux_len = length(agents);
redP   = zeros(aux_len,1);
redPhi = zeros(aux_len,1);

redP = [agents.Pm]';
redPhi = [agents.Phi]';

red_del = [];
red_ind = 1;
for j=1:aux_len
    for k=j+1:aux_len
        if (abs(redP(j)-redP(k)) < REDUNDANCY_P_MAX) & (abs(redPhi(j)-redPhi(k)) < REDUNDANCY_PHI_MAX)
            red_del(red_ind,1)=j;
            red_del(red_ind,2)=k;
            red_ind = red_ind+1;
        end
    end
end
red_del_final = zeros(length(agents),1);
for j=1:size(red_del,1)
    if agents(red_del(j,1)).age > agents(red_del(j,2)).age
        red_del_final(red_del(j,1)) = 1;
    else
        red_del_final(red_del(j,2)) = 1;
    end
end
j = 1;
while j<=length(red_del_final)
    if red_del_final(j)
        agents(j) = [];
        red_del_final(j) = [];
        KILLED_BY_REDUNDANCY = KILLED_BY_REDUNDANCY + 1;
    else
        j=j+1;
    end
end