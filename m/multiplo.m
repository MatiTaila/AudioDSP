function [n] = multiplo(Pi,Pj)


n_Cluster_Width=10;


nij=Pi/Pj;
if abs(nij-round(nij))<=n_Cluster_Width/Pj
    n=round(nij);
else
    nji=Pj/Pi;
    if abs(nji-round(nji))<=n_Cluster_Width/Pi
        n=round(nji);
    else
        n=0;
    end;
end;


end