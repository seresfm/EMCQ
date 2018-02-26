 
%Prueba para validar Throughput EMCA

clc
clear all

seed = 479;
N = 100;
NoMBS = 2;

p = 0.1;
r1 = 0.5;
r2 = 0.2;
tau_type = 9;
K = 500;


rng(seed);

A = zeros(N,N);
x = rand(N,1);
y = rand(N,1);

MBS = randsample(N, NoMBS);


for i = 1:length(A)
    for j = i+1:length(A)
        random = rand();
        dist = sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2);
        if (random <= p)
            MBS1 = length(find(MBS == i));
            MBS2 = length(find(MBS == j));
            if (MBS1== 1) || (MBS2 == 1)
                if (dist <= r1)
                    A(i,j) = 1;
                    A(j,i) = 1;
                end
            else
                if (dist <= r2)
                    A(i,j) = 1;
                    A(j,i) = 1;
                end
            end
        end
    end
end


  G = graph(A);
  plot(G)


% C?lculo de MMLR y MNE

f_e = zeros(N,N);
Num_nodos_cam = 0;
Max_Length = 0;
Node_usage = zeros(1,N);
distancias_minimas = zeros(N,N);
Minimum_Path= {};

DG = sparse(A);

for i = 1:N
    for j = 1:N
 
        if( i ~= j)
            [e, L] = graphshortestpath(DG, i, j);
            %[e L] = dijkstra(A,i,j)
            distancias_minimas(i,j) = e;
            Minimum_Path{i,j} = L;

            for ii = 1:size(L, 2)-1

            f_e(L(ii),L(ii+1)) = f_e(L(ii),L(ii+1)) + 1;
            f_e(L(ii+1),L(ii)) = f_e(L(ii+1),L(ii)) + 1;
            
            Node_usage(L(ii)) = Node_usage(L(ii)) + 1;
            if ii == (size(L, 2)-1)
            Node_usage(L(ii+1)) = Node_usage(L(ii+1)) + 1; 
            end
            
%             if ii > 1
%             Node_usage(L(ii)) = Node_usage(L(ii)) + 1; 
%             end

            end
        
           %Edges ;

           Num_nodos_cam = Num_nodos_cam + e;

           if (e > Max_Length)
               Max_Length = e;
           end
           
        end

    end
    i;
end

Num_Total_Cams = N*(N-1);

%MMLR

%MMLR = Max_Length

%MNE
Pr_e = zeros(1,N*(N-1)/2);

% for i = 1:N
%     for j =i+1:N
%        
%         Epsilon_1 = Edges(i,j) ./ Num_nodos_cam;
%         
%         Epsilon = [Epsilon  Epsilon_1];
%         
%     end
% end
k = 1;
for i = 1:N
    for j =i+1:N
       
        Pr_e(k) = f_e(i,j) ./ Num_Total_Cams;
        
        k = k+1;
        
    end
end


%Normalizaci?n de Pr
sum_f_e = sum(Pr_e);
Pr_e_norm = Pr_e ./ sum_f_e;

H_G = 0;

for i = 1:size(Pr_e_norm,2)
    
    if(Pr_e_norm(i) > 0)
        
        H_G = H_G - Pr_e_norm(i)*log2(Pr_e_norm(i));
    
    end
        
end
    
%H_G


%C?lculo de la Probabilidad de Tx

Ne = size(Pr_e_norm,2);
Hmax_e = -log2(1/Ne);
Hmax_n = -log2(1/N);
tau = zeros(N,1);
H_G_nodo = 0;
Pr_nodo = Node_usage ./ Num_Total_Cams;  %Probabilidad de uso de los nodos
sum_Pr_nodo = sum(Pr_nodo);
Pr_nodo_norm = Pr_nodo ./ sum_Pr_nodo;

for i = 1:size(Pr_nodo_norm,2)
    
    if(Pr_nodo_norm(i) > 0)
        
        %H_G_nodo = H_G_nodo - Pr_nodo(i)*log2(Pr_nodo(i));
        H_G_nodo = H_G_nodo - Pr_nodo_norm(i)*log2(Pr_nodo_norm(i));
    
    end
        
end
    
%H_G_nodo


  %Asignaci?n de la Probabilidad de Transmisi?n
    
  gamma1 = 0.7;
  gamma2 = 0.8;
  gamma3 = 0.3;
  
 if (tau_type == 1)   
    tau = exp(-(N)*(Hmax_n/H_G_nodo).*Pr_nodo_norm);   
 elseif (tau_type == 2)
     tau = exp(-(N)*(H_G_nodo/Hmax_n).*Pr_nodo_norm);
 elseif (tau_type == 3)
     tau = exp(-(Hmax_n/H_G_nodo).*Pr_nodo);
 elseif (tau_type == 4)
     tau = exp(-(H_G_nodo/Hmax_n).*Pr_nodo);
 elseif (tau_type == 5)
     tau = 0.4*ones(1,N);
 elseif (tau_type == 6)
      tau = 0.2*ones(1,N);         
 elseif (tau_type == 7)
      tau = 0.3*ones(1,N);              
 elseif (tau_type == 8)
      tau = 0.1*ones(1,N);  
 elseif (tau_type == 9)
      tau_h = 0.5*ones(1,N); 
      tau_l = 0.3*ones(1,N); 
 elseif (tau_type == 10)
      tau = 0.01*ones(1,N); 
 elseif (tau_type == 11)
      tau = 0.8*ones(1,N);
 elseif (tau_type == 12)
      tau = gamma1*exp(gamma1.*Pr_nodo_norm);
 elseif (tau_type == 13)
      tau = gamma2*exp(gamma2.*Pr_nodo_norm);
 elseif (tau_type == 14)
       tau = gamma3*exp(-gamma3.*Pr_nodo_norm.*N);
%      tau = gamma3*exp(-gamma3.*Pr_nodo)
 end

%tau

% EMCA Application

%Generaci?n de tr?fico por slots 

P_free = zeros(N,1);
P_success = zeros(N,1);
P_rx_bussy = zeros(N,1);
P_collision = zeros(N,1);
e_used_paths = zeros(N,N);
% lambda_i = zeros(N,1);

for i = 1:N
    for j = 1:N
        
        if (f_e(i,j) > 0 )
            
            e_used_paths(i,j) = 1;
            e_used_paths(j,i) = 1;
            
        end
    end
end

e_with_node = cell(N,1);
Q = cell(N,1);
Cola_h = cell(N,1);
Cola_l = cell(N,1);

% for i = 1:N
%     
%     e_with_node{i} = find(e_used_paths(i,:));
%     
% end

for i = 1:N
    nodos_e = find(e_used_paths(i,:));
    
    for k = 1:size(nodos_e,2)
        new_elements = f_e(i, nodos_e(k));
        e_with_node{i} = [e_with_node{i} repmat(nodos_e(k),1,new_elements)];       
    end
    
end

 NoSlots = 0;
 Intentos_Tx_Paqs = 0;
flag_h = zeros(N,1);
flag_l = zeros(N,1);
Pa0 = 0;
Pa1_h = 0;
Pa1_l = 0;
Pa_plus_h = 0;
Pa_plus_l = 0;
Pii_simulation = zeros(K+1,K+1);
P_states_xy = zeros(K + 1,K + 1);
Pp = 0;
m1 = 1;
m2 = 1;

   for i = 1:N    
       
     Cola_h{i} = [];
     Cola_l{i} = [];
    
   end
   
   
   %Selecci?n de un nodo promedio
   
  [Val, Ave_Node] = min(abs(Node_usage - mean(Node_usage)));
 

%tic
while Intentos_Tx_Paqs < 1000000
    Intentos_Tx_Paqs
      NoSlots = NoSlots + 1;
   for i = 1:N     
     Q{i} = [];
     flag_h(i) = 0;
     flag_l(i) = 0;
     
   end
    
    
   for i = 1:N    %Ciclo para determinar si en el slot transmite el nodo i
       
       if( tau_h(i) >= rand(1))           
           flag_h(i) = 1;   
       
       elseif (tau_l(i) >= rand(1))
           flag_l(i) = 1; 
       end   
       
       if ( i == Ave_Node ) && (isempty(Cola_h{Ave_Node}) == 1)
           
           flag_h(Ave_Node) = 0; 

       end
       
       if ( i == Ave_Node ) && (isempty(Cola_l{Ave_Node}) == 1)
           
           flag_l(Ave_Node) = 0; 

       end
       
   end
%    
%    for i = 1:N   %Ciclo para ejecutar el proceso de env?o de paquetes para todos los nodos
%        
%        if (flag(i) == 1)   
%            
%            if (size(e_with_node{i},2) > 0)
%                
%            node_j = 0;
%            while node_j < 1
% 
%            node_j = round(rand(1)*size(e_with_node{i},2));  
%            end        
%            sig_node = e_with_node{i}(node_j);
%            Q{sig_node}(end + 1) = i;
%            %sig_node
%            end
%        
%        end
%        
%    end
   
   for i = 1:N   %Ciclo para ejecutar el proceso de env?o de paquetes para todos los nodos
       
                  
               if ( i == Ave_Node)           %if para encontrar el valor de las probailidades de los estados de la cadena
               m1 = size(Cola_h{Ave_Node}, 2) + 1;
               Pii_simulation(m1,m2) = Pii_simulation(m1,m2) + 1;
               
               m2 = size(Cola_l{Ave_Node}, 2) + 1;
               Pii_simulation(m1,m2) = Pii_simulation(m1,m2) + 1;
               end
               
 %Proceso de Transmisi?n de un paquete de alta prioridad
     
       if (flag_h(i) == 1)   
           
           if (isempty(Cola_h{i}) == 1 )    % Si la cola del buffer esta vac?a, se genera un nuevo paquete a transmitir
                
               destino = 0;
               ttcuenta = 0;
               
                    while destino == 0 
                        
                        destino = round((N-1)*(rand(1)) + 1);
                        ttcuenta = ttcuenta + 1;
                        Es_high = 0;
                        for aux10 = 1:size(MBS,1)
                            aux11 = find(Minimum_Path{i, destino} == MBS(aux10),1);
                            if (isempty(aux11)==1)
                            aux11 = 0;
                            end
                        Es_high = Es_high + aux11;
                        end
                        
                        if destino == i 
                            destino = 0;
                        end
                        if Es_high == 0 
                            destino = 0;
                        end
                        if ttcuenta == 20
                            break
                        end
                    end

                   if  (i ~= Ave_Node) 

                       if (ttcuenta < 10)
%                            Minimum_Path{i, destino}
                        if(isempty(Minimum_Path{i, destino})==0)
               sig_node = Minimum_Path{i, destino}(2);
               Q{sig_node}(end + 1) = i*1000 + destino;        
               Cola_h{i} = i*1000 + destino;
                        end
                       end
                   end
                
               
           else                        %Si la cola tiene al menos un elemento, 
                                         %entonces se env?a el sig paquete que se encuentre al principio
                                         %de la cola.  
                    
               sig_node = Cola_h{i}(1);
               origen = fix(sig_node /1000);
               destino = int8((sig_node/1000 -  origen)*1000);
               
               while i == destino
                   Cola_h{i}(1) = [];
                   if (isempty(Cola_h{i}) == 0)
                   sig_node = Cola_h{i}(1);
                   origen = fix(sig_node /1000);
                   destino = (sig_node/1000 -  origen)*1000;
                   destino = int8(destino);
                   else
                       break
                   end
               end
                   if (isempty(Cola_h{i}) == 0)
                        destino = int8(destino);
                        sig = find(Minimum_Path{origen,destino}==i,1) + 1;
                        sig_node = Minimum_Path{origen,destino}(sig);
                        Q{sig_node}(end + 1) = Cola_h{i}(1);
                   end
         
           end          

       end
  
       if (flag_l(i) == 1)   
           
           if (isempty(Cola_l{i}) == 1 )    % Si la cola del buffer esta vac?a, se genera un nuevo paquete a transmitir
                
               destino = 0;
               ttcuenta = 0;
               
               if  (i ~= Ave_Node) 
                   
                    while destino == 0 
                        destino = round((N-1)*(rand(1)) + 1);
                        ttcuenta = ttcuenta + 1;
                        Es_high = 0;
                        for aux10 = 1:size(MBS,1)
                            aux11 = find(Minimum_Path{i, destino} == MBS(aux10),1);
                            if (isempty(aux11)==1)
                            aux11 = 0;
                            end
                        Es_high = Es_high + aux11;
                        end
                        
%                         if i == 90
%                             i
%                             Cola_l{i}
%                             ttcuenta
%                             destino
%                             Es_high
%                         end
                        
                        if destino == i 
                            destino = 0;
                        end
                        if Es_high >= 1 
                            destino = 0;
                        end
                        
                        if ttcuenta == 20
                            break
                        end
                    end


                       if (ttcuenta < 10)
%                            Minimum_Path{i, destino}
                if(isempty(Minimum_Path{i, destino})==0)
               sig_node = Minimum_Path{i, destino}(2);
               Q{sig_node}(end + 1) = i*1000 + destino;          
               Cola_l{i} = i*1000 + destino;
                end
                       end
                   end
                
               
           else                        %Si la cola tiene al menos un elemento, 
                                         %entonces se env?a el sig paquete que se encuentre al principio
                                         %de la cola.  
                    
               sig_node = Cola_l{i}(1);
               origen = fix(sig_node /1000);
               destino = int8((sig_node/1000 -  origen)*1000);
               
               while i == destino
                   Cola_l{i}(1) = [];
                   if (isempty(Cola_l{i}) == 0)
                   sig_node = Cola_l{i}(1);
                   origen = fix(sig_node /1000);
                   destino = (sig_node/1000 -  origen)*1000;
                   destino = int8(destino);
                   else
                       break
                   end
               end
                   if (isempty(Cola_l{i}) == 0)
                        destino = int8(destino);
                        sig = find(Minimum_Path{origen,destino}==i,1) + 1;
                        sig_node = Minimum_Path{origen,destino}(sig);
                        Q{sig_node}(end + 1) = Cola_l{i}(1);
                   end

           end          

       end
       
   
   end

%    Cola_h
%     Q

    % Ciclo para determinar si hay tx exitosas, colisiones, etc.
    
    
   for i = 1:N   
%        size(Cola{i},2)
%          lambda_i(i) = lambda_i(i) + size(Q{i},2);
   %Cola_h{90}
%     Q

        Intentos_Tx_Paqs = Intentos_Tx_Paqs + size(Q{i},2);
        aux12 = flag_h(i) + flag_l(i);
        
       if (size(Q{i},2) == 0) %&& (flag(i) == 0)

           P_free(i) = P_free(i) + 1;
           
           if (i == Ave_Node)
               
               Pa0 = Pa0 + 1;
               
           end           
           
       elseif (size(Q{i},2) == 1 && flag_l(i) == 0 && flag_h(i) == 0)           
           
           source_aux = fix(Q{i}/1000);
           sink_aux = int8((Q{i}/1000 -  source_aux)*1000);
           Source = find(Minimum_Path{source_aux, sink_aux}==i,1) - 1;
           Source = Minimum_Path{source_aux, sink_aux}(Source);
           P_success(Source) = P_success(Source) + 1;
           
           tamano_cola = size(Cola_h{i},2) + size(Cola_l{i},2);
          if (tamano_cola == K) && (i == Ave_Node)   %Condicion para suma del calculo de la probilidad de Paq perdidos
              
               Pp = Pp + 1;           
                
          end
       
                  sig_node = Q{i};
                  origen = fix(sig_node /1000);
                  destino = int8((sig_node/1000 -  origen)*1000);
                  Es_high = 0;
                  for aux10 = 1:size(MBS,1)
                      aux14 = find(Minimum_Path{origen, destino}==MBS(aux10),1);
                    if (isempty(aux14)==1)
                        aux14 = 0;
                    end
                  Es_high = Es_high + aux14;
                  end
                        
                        
           if (tamano_cola < K) %&& (i == Ave_Node)  %Condicion que limita el tama?o del buffer de un nodo a K
                
           if(Es_high > 0)
               Cola_h{i}(end + 1) = Q{i};
               Cola_h{Source}(1) = [];
               else
               Cola_l{i}(end + 1) = Q{i};
               Cola_l{Source}(1) = [];
           end
          
           
           if (i == Ave_Node && Es_high > 0)
               
               Pa1_h = Pa1_h + 1;
           
           elseif(i == Ave_Node && Es_high == 0)
               Pa1_l = Pa1_l + 1;
           end
           
           end
           
       elseif (size(Q{i},2) > 1 && flag_h(i) == 0 &&flag_l(i) == 0)
           j = 1;
           while j <= size(Q{i},2)
                source_aux = fix(Q{i}(j)/1000);
                sink_aux = int8((Q{i}(j)/1000 -  source_aux)*1000);
                Source = find(Minimum_Path{source_aux, sink_aux}==i,1) - 1;
                Source = Minimum_Path{source_aux, sink_aux}(Source);
                
               P_collision(Source) = P_collision(Source) + 1;
               j = j + 1;
           end
           
           for j = 1:size(Q{i},2)
               
                  sig_node = Q{i}(j);
                  origen = fix(sig_node /1000);
                  destino = int8((sig_node/1000 -  origen)*1000);
                  Es_high = 0;              
                  for aux10 = 1:size(MBS,1)
                      aux14 = find(Minimum_Path{origen, destino}==MBS(aux10),1);
                    if (isempty(aux14)==1)
                        aux14 = 0;
                    end
                    Es_high = Es_high + aux14;
                  end
           end
           Es_low = size(Q{i},2) - Es_high;
                  
           if (i == Ave_Node && Es_high > 1)
               
               Pa_plus_h = Pa_plus_h + 1;
               
           elseif (i == Ave_Node && Es_low > 1 )
               
               Pa_plus_l = Pa_plus_l + 1;
               
           end
           
               
       elseif (size(Q{i},2) >= 1 && aux12 >= 1 )
           j = 1;
           while j <= size(Q{i},2)
                source_aux = fix(Q{i}(j)/1000);
                sink_aux = int8((Q{i}(j)/1000 -  source_aux)*1000);
                Source = find(Minimum_Path{source_aux, sink_aux}==i,1) - 1;
                Source = Minimum_Path{source_aux, sink_aux}(Source);
               P_rx_bussy(Source) = P_rx_bussy(Source) + 1;
               j = j + 1;
           end

           if size(Q{i},2) == 1
                  sig_node = Q{i};
                  origen = fix(sig_node /1000);
                  destino = int8((sig_node/1000 -  origen)*1000);
                  Es_high = 0;
                  for aux10 = 1:size(MBS,1)
                      aux14 = find(Minimum_Path{origen, destino}==MBS(aux10),1);
                    if (isempty(aux14)==1)
                        aux14 = 0;
                    end
                  Es_high = Es_high + aux14;
                  end
           end
           
           if (size(Q{i},2) == 1) && (i == Ave_Node) && (Es_high > 0)
               
               Pa1_h = Pa1_h + 1;
               
           elseif (size(Q{i},2) == 1) && (i == Ave_Node) && (Es_high ==0)
               
               Pa1_l = Pa1_l + 1;
               
           elseif (size(Q{i},2) > 1) && (i == Ave_Node) && (Es_high > 0)
               
               Pa_plus_h = Pa_plus_h+ 1;
               
           elseif (size(Q{i},2) > 1) && (i == Ave_Node) && (Es_high ==0)
               
               Pa_plus_l = Pa_plus_l + 1;
               
           end
       end
      

   end
 
%                l = size(Cola{Ave_Node}, 2) + 1;
%                
%                P_states_xy (m, l) = P_states_xy (m, l) + 1;
               

    
end

%toc


% Pk_k = 0;
% Pk_k1 = 0;
% Pk_k_1 = 0;
% 
% for i = 1:K+1
%     for j = 1:K+1
%         
%         if (i == j)          
%             Pk_k = Pk_k + P_states_xy(i, j);
%         elseif(i > j)
%             Pk_k1 = Pk_k1 + P_states_xy(i, j);
%         elseif (i < j)
%             Pk_k_1 = Pk_k_1 + P_states_xy(i, j);
%         end
%     end
% end

% Pk_k = Pk_k / NoSlots;
% Pk_k1 = Pk_k1 / NoSlots;
% Pk_k_1 = Pk_k_1 / NoSlots;

sum_Pii = sum(sum(Pii_simulation));
Pii_simulation = Pii_simulation / sum_Pii;

Psuccess = sum(P_success);
Pcollision = sum(P_collision);
Pfree = sum(P_free);
Prxbussy = sum(P_rx_bussy);

Ptx = Psuccess + Pcollision + Pfree + Prxbussy;

Psuccess = Psuccess/Ptx;
Pcollision = Pcollision/Ptx;
Pfree = Pfree/Ptx;
Prxbussy = Prxbussy/Ptx;

Throughput = sum(P_success)/NoSlots;

La = 0;
Lb = 0;
L_T = 0;

for i= 1: K
    for j = 1: K+1-i
    
    La = (i-1)*Pii_simulation(i,j) + La;
    
    end
end

for i= 1: K
    for j = 1: K+1-i
    
    Lb = (j-1)*Pii_simulation(i,j) + Lb;
    
    end
end

for i= 1: K
    for j = 1: K+1-i
    
    L_T = (i + j)*Pii_simulation(i,j) + L_T;
    
    end
end

Pa1_h = Pa1_h / NoSlots;
Pa_plus_h = Pa_plus_h / NoSlots;
Pp = Pp / NoSlots;

Da = La / (Pa1_h + Pa_plus_h);

Pa1_l = Pa1_l / NoSlots;
Pa_plus_l = Pa_plus_l / NoSlots;
Pp = Pp / NoSlots;

Db = Lb / (Pa1_l + Pa_plus_l);

D_T = L_T / (Pa1_h + Pa1_l);
