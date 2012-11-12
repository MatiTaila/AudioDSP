function [rn] = r(Pi,Pj,n_hop,fs)

hop    = n_hop/fs; 
n_Cluster_Width=fix((25*10^(-3))/hop);


nij=Pi/Pj;
if abs(nij-round(nij))<=n_Cluster_Width/Pj
    n=round(nij);
else nji=Pj/Pi;
    if abs(nji-round(nji))<=n_Cluster_Width/Pi
        n=round(nji);
    else n=0;
    end;
end;

if (1<=n)&(n<=8)
    if (1<=n)&(n<=4)
        rn=6-n;
    else rn=1;
    end;
else rn=0;
end;


