function [agents] = tracking(agents,MaxTabSF,n_Pmax,n_Pmin,hop,MAX_AGENTS,MAX_OUTER,opt,REDUNDANCY_P_MAX,REDUNDANCY_PHI_MAX)

id = length(agents)+1;

KILLED_BY_LOSS = 0;
KILLED_BY_REDUNDANCY = 0;
KILLED_BY_OBSOLENCE = 0;
KILLED_BY_REPLACEMENT = 0;
KILLED_BY_OUTER_REGION = 0;

for i=1:size(MaxTabSF,1)-1
%     keyboard
    delete = zeros(length(agents),1);
    inner_count = 0;
    outer_count = 0;
    waiting_agents = [];
    ind = 1;
    for j=1:length(agents)
        Tout_R = 0.4*agents(j).Pm(end);
        Tout_L = 0.2*agents(j).Pm(end);
        Tin    = round(46.4e-3/hop);
        m      = MaxTabSF(i,1);
        if agents(j).Phi(end)~=0
            bp = agents(j).Phi(end)+agents(j).Pm;
        else
            bp = agents(j).pre;
        end
        error  = m-bp;
        
        while bp+Tout_R < m
            bp=bp+agents(j).Pm(end);
            error = m-bp;
        end
        
        if(abs(MaxTabSF(i+1,1)-bp)>=abs(error)) %& (abs(MaxTabSF(i-1,1)-bp)>=abs(error))
            if ( error > Tout_R) || ( error < -Tout_L )
                delete(j) = 1;
                agents(j).outer_loss = agents(j).outer_loss + 1;
            else % cae dentro de alguno de los intevalos
                if abs(error)<Tin % inner region
                    inner_count = inner_count+1;
                    agents(j).loss = 0;
                    agents(j).outer_loss = 0;
                    agents(j).Pm  = agents(j).Pm(end)+0.25*error;
                    % limits for P
                    if agents(j).Pm(end)>n_Pmax;
                        agents(j).Pm(end) = n_Pmax;
                        if opt.log>=1
                            disp('inner: saturacion periodo por arriba')
                        end
                    elseif agents(j).Pm(end)<n_Pmin
                        agents(j).Pm(end)=n_Pmin;
                        if opt.log>=1
                            disp('inner: saturacion periodo por abajo')
                        end
                    end
                    if agents(j).Phi(end)~=0
                        agents(j).Phi = [agents(j).Phi agents(j).Phi(end)+agents(j).Pm(end)];
                    else
                        agents(j).Phi = bp+0.25*error;
                    end
                    delta_s = (1-abs(error)/Tout_R)*MaxTabSF(i,2)*agents(j).Pm(end)/n_Pmax;
                else % outer region
                    % creo hijos
                    outer_count = outer_count+1;
                    agents(j).loss = agents(j).loss + 1;
                    agents(j).outer_loss = 0;
                    agents(j).Pm(end) = agents(j).Pm(end);
                    % limits for P
                    if agents(j).Pm(end)>n_Pmax;
                        agents(j).Pm(end) = n_Pmax;
                        if opt.log>=1
                            disp('outer: saturacion periodo por arriba')
                        end
                    elseif agents(j).Pm(end)<n_Pmin
                        agents(j).Pm(end)=n_Pmin;
                        if opt.log>=1
                            disp('outer: saturacion periodo por abajo')
                        end
                    end
                    if agents(j).Phi(end)~=0
                        agents(j).Phi = [agents(j).Phi agents(j).Phi(end)+agents(j).Pm(end)];
                    else
                        agents(j).Phi = bp+0.25*error;
                    end
                    delta_s = -(abs(error)/Tout_R)*MaxTabSF(i,2)*agents(j).Pm(end)/n_Pmax;
                    P_hijos = {0,error,0.5*error};
                    Phi_hijos = {error,error,0.5*error};
                    for k=1:3
                        waiting_agents(ind).pre = 0;
                        waiting_agents(ind).Pm = agents(j).Pm(end)+P_hijos{k};
                        % limits for P
                        if waiting_agents(ind).Pm(end)>n_Pmax;
                            waiting_agents(ind).Pm(end) = n_Pmax;
                            if opt.log>=1
                                disp('hijo: saturacion periodo por arriba')
                            end
                        elseif waiting_agents(ind).Pm(end)<n_Pmin
                            waiting_agents(ind).Pm(end)=n_Pmin;
                            if opt.log>=1
                                disp('hijo: saturacion periodo por abajo')
                            end
                        end
                        if agents(j).Phi(end)~=0
                            waiting_agents(ind).Phi = [agents(j).Phi agents(j).Phi(end)+Phi_hijos{k}+waiting_agents(ind).Pm(end)];
                        else
                            waiting_agents(ind).Phi = bp+Phi_hijos{k}+waiting_agents(ind).Pm(end);
                        end
                        waiting_agents(ind).Sraw = 0;
                        waiting_agents(ind).Srel = 0;
                        waiting_agents(ind).S = 0.9*agents(j).S(end);
                        waiting_agents(ind).age = 0;
                        waiting_agents(ind).loss = 0;
                        waiting_agents(ind).outer_loss = 0;
                        waiting_agents(ind).wins = 0;
                        ind = ind+1;
                    end
                end
                agents(j).S = agents(j).S+delta_s;
            end    
        else
            agents(j).S = agents(j).S*.99;
        end
        agents(j).age=agents(j).age+1;
    end
    
    if opt.log >= 2
        fprintf('pico numero: %d\n',i)
        fprintf('cantidad de agentes en el inner: %d\n',inner_count)
        fprintf('cantidad de agentes en el outer: %d\n',outer_count)
        fprintf('---------------------------------------\n')
    end
    
    % delete agents / capaz q es mejor no borrarlos sino guardarlos para
    % conservar toda la historia e ir viendo cuando mueren
    j = 1;
    while j<=length(agents)
        if delete(j) & (agents(j).outer_loss>MAX_OUTER) % MAX_OUTE R is used both when bp is in outer region or outside inner and outer
            agents(j) = [];
            delete(j) = [];
            KILLED_BY_OUTER_REGION = KILLED_BY_OUTER_REGION + 1;
        else
            j=j+1;
        end
    end
    
    % redundancy
    if exist('REDUNDANCY_P_MAX','var')
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
    end
    
    % loss
    j=1;
    while j<=length(agents)
        if agents(j).loss > MAX_OUTER
            agents(j) = [];
            KILLED_BY_LOSS = KILLED_BY_LOSS + 1;
        else
            j=j+1;
        end
    end
    
    % obsolescence
    aux = zeros(length(agents),1);
    for j=1:length(agents)
        aux(j)=agents(j).S(end);
    end
    [maxS,index] = max(aux);
    j=1;
    while j<=length(agents)
        if abs(agents(j).S(end)-maxS) > 0.8*maxS
            agents(j) = [];
            KILLED_BY_OBSOLENCE = KILLED_BY_OBSOLENCE + 1;
        else
            j=j+1;
        end
    end
    
    % replacement
    while length(agents)>MAX_AGENTS
        aux = zeros(length(agents),1);
        for j=1:length(agents)
            aux(j)=agents(j).S(end);
        end
        [a,index] = min(aux);
        agents(index)=[];
        KILLED_BY_REPLACEMENT = KILLED_BY_REPLACEMENT + 1;
    end
    
    % add waiting agents
    L = length(agents);
    for j=1:length(waiting_agents)
        waiting_agents(j).id = id;
        id = id+1;
        agents(L+j)=waiting_agents(j);
    end
    
    if opt.log >= 1
        if length(agents) == 0
            fprintf('mierda: %d\n=======================================\n',i)
        end
    end
    
    % Alternative Referee
    Svec = [agents.S]';
    [unUsed,Sindex] = max(Svec);
    agents(Sindex).wins = agents(Sindex).wins + 1;
    
end

if opt.log>=1
    fprintf('-----------------------------------------------\n')
    fprintf('KILLED\n')
    fprintf('-----------------------------------------------\n')
    fprintf('KILLED_BY_OUTER_REGION: %d\n',KILLED_BY_OUTER_REGION)
    fprintf('KILLED_BY_REDUNDANCY: %d\n',KILLED_BY_REDUNDANCY)
    fprintf('KILLED_BY_OBSOLENCE: %d\n',KILLED_BY_OBSOLENCE)
    fprintf('KILLED_BY_LOSS: %d\n',KILLED_BY_LOSS)
    fprintf('KILLED_BY_REPLACEMENT: %d\n',KILLED_BY_REPLACEMENT)
    fprintf('===============================================\n')
end