function [PF, FDev]= pdfFitForce_2motors(k1,k_1,k2,k_2,k3,f0,N,Data,DataDev,setPDF)

    %%%CALCOLO PF da ME+G_f o da VK+G_f 
% 
% parSet= 'SYNTH-A1A2-N20-f3-T25';
% fileName= '../datiGillespie/20240515/SYNTH-A1A2-N20-f3-T25/datiIstoGillForce_SYNTH-A1A2-N20-f3-T25.dat';
% Data= load(fileName);
% fileNameDevT= sprintf('../datiGillespie/%s/%s/datiGill_forceDev-av_%s_it800.dat',expNum,parSetT,parSetT);
% 
%         ForceDev= load(fileNameDev);
%         %%Y= 0;%0.5*randn(length(ForceDev));
%         %plot(DataDev(:,1), DataDev(:,2))
%         DataDev= [ForceDev(:,1) ForceDev(:,2)];%0.05*rand(1)];
%         %plot(DataDev(:,1), DataDev(:,2)),shg
%         %DataDev= DataDev(1:300,:);
%         DataDev= DataDev(1:4000,:);
% 
% % %Costanti di reazione    
% [k1, k_1, k2, k_2, k3, f0, N]= readReactionConst(parSet);
% setPDF= 'ME';
% %setPDF= 'VK';

addpath 'C:\Users\Vali\Desktop\MATLAB\PhD\HMMnanomachine\StochSimulations'

switch setPDF
    case 'ME'
        %%%%Calcolo della P(n1,n2)
        [Q, q1, q2]= masterEq(N, k1,k_1,k2,k_2,k3);

        m= 1:((N+1)*(N+2))/2;
        V= null(Q-diag(sum(Q,1)));
        V= V/sum(V);

        ME= [q1 q2 m' V];

        %%%%Calcolo della pdf della forza F
        x= Data(:,1);
        TRESH= 10^-4;

        PF= zeros(length(x),1);
        for l= 2:((N+1)*(N+2))/2
            if ME(l,4)< TRESH
                continue
            else
                qq1= ME(l,1)-1;
                qq2= ME(l,2)-1;

                MU= qq2*((f0*11)/20);
                SIGMA2= (f0^2)*((qq1/3)+ qq2*(27/400));

                PF_temp= (1/(sqrt(2*pi*SIGMA2)))*exp(-(x-MU).^2/(2*SIGMA2));
                %[qq1 qq2 ME(l,4)]
                %   figure(1)
                % plot(Data(:,1), PF_temp*ME(l,4)), hold on
                % sum(PF_temp)*(Data(3,1)-Data(2,1))

                %PdfSum= PdfSum+ PF_temp*M(l,4);
                PF= PF+ PF_temp*ME(l,4);
            end
        end
        %     figure(2)
        % bar(Data(:,1),Data(:,2),'FaceColor','0.00,0.45,0.74','EdgeColor','k','LineWidth',0.5,'BarWidth',1),hold on
        % plot(Data(:,1),PF,'Color','0.85,0.33,0.10','LineWidth',1),shg
        % 
    case 'VK'
            %%%%Calcolo della P(n1,n2)
        G= (k_1*(k_2+k3)+k2*k3)/(k2+k_2+k3);
        % punti fissi
        yStar= (k1*(k_2+k3))/((k1+G)*(k2+k_2+k3));
        zStar= (k1*k2)/((k1+G)*(k2+k_2+k3));

        M= zeros(3,3);  %matrice dei coefficienti
        b= zeros(3,1);  %vettore colonna dei termini noti
        % matrice jacobiana valutata all'equilibrio
        J= zeros(2,2);
        J(1,1)= -(k1+ k_1+ k2);
        J(1,2)= -k1+ k_2;
        J(2,1)= k2;
        J(2,2)= -(k_2+ k3);

        % momenti secondi della distro delle fluttuazioni
        M(1,1)= 2*J(1,1);
        M(1,2)= 0;
        M(1,3)= 2*J(1,2);
        M(2,1)= 0;
        M(2,2)= 2*J(2,2);
        M(2,3)= 2*J(2,1);
        M(3,1)= J(2,1);
        M(3,2)= J(1,2);
        M(3,3)= J(1,1)+J(2,2);

        b(1)= k1*(1- yStar- zStar)+ (k_1+ k2)*yStar+ k_2*zStar; % B11
        b(2)= k2*yStar+ (k_2+ k3)*zStar;                        % B22
        b(3)= -k2*yStar- k_2*zStar;                             % B12

        zeta= -inv(M)*b;
        xi2= zeta(1);
        eta2= zeta(2);
        xieta= zeta(3);

        rho= xieta/(sqrt(xi2*eta2));
        Sigma= xi2*eta2*(1-rho^2);
       
        q= 0.05:0.5:N;    
        [qq1, qq2]= meshgrid(q);

        Pn1n2= 1/(2*pi*sqrt(N^2*Sigma))*exp(-(1/(2*(1-rho^2)*N))*(((qq1-N*yStar).^2)/xi2+ ((qq2-N*zStar).^2)/eta2- 2*rho.*(qq1-N*yStar).*(qq2-N*zStar)./(sqrt(xi2*eta2))));

% figure()
% surf(qq1,qq2,Pn1n2)


        %maxf= max(Data(:,1)); %20;
        %minf= min(Data(:,1)); %-10;
        %df= Data(3,1)-Data(2,1);
        %effe= Data(:,1);
        index=1;
        %PF= zeros(length(minf:df:maxf)+1,1);
        PF= zeros(length(Data),1);
   
        %for effe= minf:df:maxf
        for i= 1:length(Data)
            effe= Data(i,1);

            MU= qq2.*((f0*11)/20);
            SIGMA2= (f0^2)*((qq1./3)+ qq2.*(27/400));
            PF_temp= (1./(sqrt(2*pi*SIGMA2))).*exp(-(effe-MU).^2./(2*SIGMA2));

            %figure()
            %surf(qq1,qq2,Pn1n2.*PF_temp)

            PF(index)= trapz(q, trapz(q, Pn1n2.*PF_temp,2));
            index= index+1;

            %figure()
            %surf(PF_temp)
            %view(2);
            %shading interp;
            %figure(50)
            %plot(effe(e), Peffe(e), '*'), hold on
        end
        % figure(50)
        % bar(Data(:,1),Data(:,2),'FaceColor','0.00,0.45,0.74','EdgeColor','k','LineWidth',0.5,'BarWidth',1),hold on
        % plot(Data(:,1),Data(:,2)), hold on
        % plot(Data(:,1), PF, '*'), hold on
end

    %%%%CALCOLO FtotDev
G= (k_1*(k_2+k3)+k2*k3)/(k2+k_2+k3);    
yStar= (k1*(k_2+k3))/((k1+G)*(k2+k_2+k3));
zStar= (k1*k2)/((k1+G)*(k2+k_2+k3));
%     %Cond iniziali
% concIn= [0 0];
% xD= 1- concIn(1)- concIn(2);
%b= zeros(2,1);
%b(1)= k1;
%b(2)= 0;

A= zeros(2);
A(1,1)=-(k1+k_1+k2);
A(1,2)= -(k1-k_2);
A(2,1)= k2;
A(2,2)= -(k_2+k3);
[V,D]= eig(A);
xStar= [yStar 
    zStar];
%xStar= -inv(A)*b;
% DET_V=V(1,1)*V(2,2)-V(1,2)*V(2,1)
c= -inv(V)*xStar;
%c1= (-xStar(1)*V(2,2)+ xStar(2)*V(1,2))/(det(V));
%c2= (xStar(1)*V(2,1)- xStar(2)*V(1,1))/(det(V));


time= DataDev(:,1)';
% tMax= 0.5;
% dt= 0.001;
% time= 0:dt:tMax;

%x= [];
x= xStar+ exp(D(1,1).*time)*c(1).*V(:,1)+ exp(D(2,2).*time)*c(2).*V(:,2);
%y= xStar(1)+ exp(D(1,1)*time)*c1*V(1,1)+ exp(D(2,2)*time)*c2*V(1,2);
%z= xStar(2)+ exp(D(1,1)*time)*c1*V(2,1)+ exp(D(2,2)*time)*c2*V(2,2);
x= real(x);
FDev= x(2,:)*(11/20)*f0*N;
FDev= FDev';

%     figure(1)
% plot(time, x(1,:), 'k'), hold on
% plot(time, x(2,:), 'r'), hold on
% %plot(time, y, '--b'), hold on
% %plot(time, z, '--y'), hold on
% xlabel('$t\ (s)$','fontsize',16)
% ylabel('$y(t), \ z(t)$','fontsize',16)
% legend('$y(t)$','$z(t)$')
% 
% 
%     figure(3)
% plot(time, FDev, '--', 'Color', '#a2142f', 'Linewidth', 1.0), hold on
% xlabel('$t\ (s)$','fontsize',16)
% ylabel('$F\ (pN)$','fontsize',16)