% gillespieAlgorithm

function [nD, nA1, nA2, F, tempo, atpTot]= gillespieAlgorithm(N, concIn, tmax, k1, k_1, k2, k_2, k3, f0)

   % Parametri simulazione
t= 0;
tempo= [];%zeros(1,tmax);
prob= 0;

    % Concentrazione iniziale
ConcA1_0= concIn(1);        
ConcA2_0= concIn(2);
%ConcD_0= 1-ConcA1_0-ConcA2_0;

    % Numero di elementi delle specie D,A1,A2
A1= floor(ConcA1_0*N);
A2= floor(ConcA2_0*N);
D= N-A1-A2; %floor(ConcD_0*N); 

    % Numero di elementi delle specie in ogni posizione al tempo 1
nD= [];%zeros(1,tmax);
nA1= [];%zeros(1,tmax);
nA2= [];%zeros(1,tmax);
   
nD(1,1)= D;
nA1(1,1)= A1;
nA2(1,1)= A2;
     
%F1= -F0+2*F0.*rand(1,tmax);  
%F2= (F0/10)+ F0*(0.9).*rand(1,tmax);  
F= []; %zeros(tmax,3);

atp= 0; %zeros(1,tmax);
atpTot=[];
TAU=[];

for i= 2:tmax

        % Tassi di transizione 
    a1= k1*D/N; %(k1/Omega)*(1- A1/N- A2/N).*ones(Omega,1); %(k1/Omega)*(1- (sum(A1+A2)/N))*ones(Omega,1);
    a2= k_1*A1/N;
    a3= k2*A1/N;
    a4= k_2*A2/N;
    a5= k3*A2/N;
    
    
    a0= a1+a2+a3+a4+a5;
    if a0<10^-17
        t
        break
    end
        
        % Algoritmo di Gillespie
    r1= rand(1,1);
    r2= rand(1,1);
    
    tau= -1/a0*log(r1); 
    r2= a0*r2;
    ind= 1;
    
    while(ind),
            prob= a1;
            if(r2< prob),
                D= D-1;
%                  if D<0
%                      t
%                      break;
%                  end
                A1= A1+1;
                ind= 0;
                
                break;
            end  
            prob= prob+ a2;
            if(r2< prob),
                A1= A1-1;
%                 if A1<0, break;
%                 end
                D= D+1;
                ind= 0;
                break;
            end   
            prob= prob+ a3;
            if(r2< prob),
                A1= A1-1;
%                 if A1<0, break;
%                 end
                A2= A2+1;
                ind= 0;
                break;
            end
            prob= prob+ a4;
            if(r2< prob),
                A2= A2-1;
%                 if A2<0, break;
%                 end
                A1= A1+1;
                ind= 0;
                break;
            end
            prob= prob+ a5;
            if(r2< prob),
                A2= A2-1;
%                 if A2<0, break;
%                 end
                D= D+1;
                ind= 0;
                atp =atp+1;
                break;
            end
    end
            
%Si aggiorna il tempo di Gillespie
    t= t+tau;
    tempo(1,i)= t;
    
    TAU=[TAU
        tau];
  
%Si aggiornano le specie    
	nD(1,i)= D;
    nA1(1,i)= A1;
    nA2(1,i)= A2;
    
%Calcolo delle forze
    if nA2(1,i)>0
        F2=(f0/10)+ f0*(0.9).*rand(nA2(1,i),1);
    else 
        F2=0;
    end
    if nA1(1,i)>0
        F1= -f0+2*f0.*rand(nA1(1,i),1);
    else
        F1=0;
    end
    
    F(i,:)=[sum(F1) sum(F2) sum(F1)+sum(F2)];
    
%Calcolo dell'ATPasi
    atpTot(i)= atp;  
end