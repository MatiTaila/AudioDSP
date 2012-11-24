function agents = tracking(agents,MaxTabSF,n_Pmax,hop,MAX_AGENTS,MAX_OUTER,opt,REDUNDANCY_P_MAX,REDUNDANCY_PHI_MAX)

KILLED_BY_LOSS = 0;
KILLED_BY_REDUNDANCY = 0;
KILLED_BY_OBSOLENCE = 0;
KILLED_BY_REPLACEMENT = 0;

for i=1:size(MaxTabSF,1)-1
    
%     if i == 45, keyboard, end
    
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
        bp     = agents(j).Phi(end);
        error  = m-bp;
        
        while bp+Tout_R < m
            bp=bp+agents(j).Pm(end);
        end
        
        if(abs(MaxTabSF(i+1,1)-bp)>=abs(error))
            
            if ( error > Tout_R) || ( error < -Tout_L )
                delete(j) = 1;
            else % cae dentro de alguno de los intevalos
                if abs(error)<Tin % inner region
                    inner_count = inner_count+1;
                    agents(j).loss = 0;
                    delta_s = (1-abs(error)/Tout_R)*MaxTabSF(i,2)*agents(j).Pm(end)/n_Pmax;
                    agents(j).Pm  = agents(j).Pm(end)+0.25*error;
%                     agents(j).Phi = agents(j).Phi(end)+agents(j).Pm(end);
%                     agents(j).Pm  = [agents(j).Pm agents(j).Pm(end)+0.25*error];
                    agents(j).Phi = [agents(j).Phi agents(j).Phi(end)+agents(j).Pm(end)];

                else % outer region
                    % creo hijos
                    outer_count = outer_count+1;
                    agents(j).loss = agents(j).loss + 1;
                    delta_s = -(abs(error)/Tout_R)*MaxTabSF(i,2)*agents(j).Pm(end)/n_Pmax;
                    
                    agents(j).Pm(end)  = agents(j).Pm(end);
%                     agents(j).Phi = agents(j).Phi(end)+agents(j).Pm(end);
%                     agents(j).Pm  = [agents(j).Pm agents(j).Pm(end)];
                    agents(j).Phi = [agents(j).Phi agents(j).Phi(end)+agents(j).Pm(end)];

                    P_hijos = {0,error,0.5*error};
                    Phi_hijos = {error,error,0.5*error};
                    for k=1:3
                        waiting_agents(ind).Pm = agents(j).Pm(end)+P_hijos{k};
%                         waiting_agents(ind).Phi = agents(j).Phi(end)+Phi_hijos{k}+waiting_agents(ind).Pm(end);
%                         waiting_agents(ind).Pm = [agents(j).Pm agents(j).Pm(end)+P_hijos{k}];
                        waiting_agents(ind).Phi = [agents(j).Phi agents(j).Phi(end)+Phi_hijos{k}+waiting_agents(ind).Pm(end)];
                        
                        waiting_agents(ind).Sraw = 0;
                        waiting_agents(ind).Srel = 0;
%                         waiting_agents(ind).S = [agents(j).S 0.9*agents(j).S(end)];
                        waiting_agents(ind).S = 0.9*agents(j).S(end);
                        waiting_agents(ind).age = 0;
                        waiting_agents(ind).loss = 0;
                        
                        ind = ind+1;
                    end
                end
%                 agents(j).S = [agents(j).S agents(j).S+delta_s];
                agents(j).S = agents(j).S+delta_s;
            end    
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
        if delete(j)
            agents(j) = [];
            delete(j) = [];
        else
            j=j+1;
        end
    end
    
    % redundancy
    if exist('REDUNDANCY_P_MAX','var')
        aux_len = length(agents);
        redP   = zeros(aux_len,1);
        redPhi = zeros(aux_len,1);
        for j=1:aux_len
            redP(j)   = agents(j).Pm(end);
            redPhi(j) = agents(j).Phi(end);
            %         if length(agents(j).Phi)>1
            %             redPhi(j) = agents(j).Phi(end)-agents(j).Phi(end-1)-agents(j).Pm(end);
            %         else
            %             redPhi(j) = agents(j).Phi(end);
            %         end
        end
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
        red_del_final = zeros(size(red_del,1),1);
        for j=1:size(red_del,1)
            if agents(red_del(j,1)).age > agents(red_del(j,2)).age
                red_del_final(red_del(j,2)) = 1;
            else
                red_del_final(red_del(j,1)) = 1;
            end
        end
        j = 1;
        while j<=length(red_del_final)
            if red_del_final(j)
                agents(red_del_final(j)) = [];
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
        agents(L+j)=waiting_agents(j);
    end
    
%     % Debug
%     P = zeros(length(agents),1);
%     S = zeros(length(agents),1);
%     for j=1:length(agents)
%         P(j)=agents(j).Pm(end);
%         S(j)=agents(j).S(end);
%     end
%     [a,b]=max(S);
%     periodos(i) = agents(b).Pm(end);
%     fases(i) = agents(b).Phi(end);
    
    if opt.log >= 1
        if length(agents) == 0
            fprintf('mierda: %d\n=======================================\n',i)
        end
    end
    
end

fprintf('-----------------------------------------------\n')
fprintf('KILLED\n')
fprintf('-----------------------------------------------\n')
fprintf('KILLED_BY_REDUNDANCY: %d\n',KILLED_BY_REDUNDANCY)
fprintf('KILLED_BY_OBSOLENCE: %d\n',KILLED_BY_OBSOLENCE)
fprintf('KILLED_BY_LOSS: %d\n',KILLED_BY_LOSS)
fprintf('KILLED_BY_REPLACEMENT: %d\n',KILLED_BY_REPLACEMENT)
fprintf('===============================================\n')