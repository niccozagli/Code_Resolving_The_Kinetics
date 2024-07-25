% SIMULAZIONI GILLESPIE NS

addpath('C:\Users\Vali\Desktop\MATLAB\PhD');
setFigure;

saveResults= 0;         %SALVARE I RISULTATI? 
     
parSet= 'SYNTH-A1A2-N20-f3-T15-prova';
T= parSet(7:9);
tmax= 10^5;

    %Costanti di reazione e forza isometrica per motore   
[k1, k_1, k2, k_2, k3, f0, N]= readReactionConst(parSet);

storeXiEta= [];
storeZeta= [];

G= (k_1*(k_2+k3)+k2*k3)/(k2+k_2+k3);
dutyRatio= k1/(k1+G);

% campo medio analitico
yStar= (k1*(k_2+k3))/((k1+G)*(k2+k_2+k3));
zStar= (k1*k2)/((k1+G)*(k2+k_2+k3));

% condizione iniziale
%concIn= [yStar, zStar];
concIn= [0, 0];

[nD, nA1, nA2, F, tempo, atpTot]= gillespieAlgorithm(N, concIn, tmax, k1, k_1, k2, k_2, k3, f0);

Force= [tempo'/N F];    %F= [F1 F2 Ftot]
Popul= [tempo'/N nD' nA1' nA2'];

    % istogrammi concentrazioni
P1=[];
P2=[];
PD=[];
for i=0:N
    l1= find(Popul(:,3)==i);
    P1(i+1)= length(l1);
    l2= find(Popul(:,4)==i);
    P2(i+1)= length(l2);
    lD= find(Popul(:,2)==i);
    PD(i+1)= length(lD);
end
%     sum(P1);
%     length(Popul);
dx= 1/N;
xConc= 0:dx:1;
P1= P1/(sum(P1)*dx);
P2= P2/(sum(P2)*dx);
PD= PD/(sum(PD)*dx);
MeanP1= sum(xConc.*P1)/sum(P1)
MeanP2= sum(xConc.*P2)/sum(P2)
MeanP0= sum(xConc.*PD)/sum(PD)

    % istogrammi forza
    
tIso= floor(length(Force)*2/3);
 
[Yf Xf]=hist(Force(tIso:end,4),30);
Yf= Yf/sum(Yf)/(Xf(2)-Xf(1));
istoF=[Xf' Yf'];
    
[Yf1 Xf1]=hist(Force(tIso:end,2),100);
Yf1= Yf1/sum(Yf1)/(Xf1(2)-Xf1(1));
istoF1=[Xf1' Yf1'];

deltaX= max(Force(:,3))/35;
[Yf2 Xf2]=hist(Force(tIso:end,3), 0:deltaX:max(Force(:,3)));
Yf2= Yf2/sum(Yf2)/(Xf2(4)-Xf2(3));
istoF2=[Xf2' Yf2'];
istoF2= istoF2(2:end,:);

MeanFtot= sum(istoF(:,1).*istoF(:,2))/sum(istoF(:,2))
MeanF1= sum(istoF1(:,1).*istoF1(:,2))/sum(istoF1(:,2))
MeanF2= sum(istoF2(:,1).*istoF2(:,2))/sum(istoF2(:,2))

        % duty ratio e ATPasi
dutyratio= mean((nA1/N)+ (nA2/N));
dutyratio_EFF= mean(nA2/N);
%nTot= D+sum(A1)+sum(A2);
ATPasi= atpTot(end)/tempo(end);
ATPtheo_EFF= k3*dutyratio_EFF;
ATPtheo= (k1*k2*k3)/((k1+G)*(k2+k_2+k3));
% fluttuazioni (Van Kampen)
xi= ((nA1/N)- yStar)*sqrt(N);
eta= ((nA2/N)- zStar)*sqrt(N);
xieta= xi.*eta;

XI2= mean((xi).*(xi));
errXI2= std((xi).*(xi));
ETA2= mean((eta).*(eta));
errETA2= std((eta).*(eta));
XIETA= mean((xi).*(eta));
errXIETA= std((xi).*(eta));

storeXiEta= [storeXiEta
    k2*ones(length(tempo'),1) tempo'/N xi' eta'];
storeZeta= [storeZeta
    k2 XI2 errXI2 ETA2 errETA2 XIETA errXIETA];

    % distribuzione di probabilità delle fluttuazioni Pi(xi, eta)
 
    % FIGURE

if(0)   % FORZA F2 E Ftot IN FUNZIONE DEL TEMPO
    figure(1)
%plot(tempo, F(:,3),'o-'), hold on
plot(storeGillForce(:,1)/N, storeGillForce(:,3),'+-'),hold on
plot(storeGillForce(:,1)/N, storeGillForce(:,2),'-'),hold on
xlabel('$t(s)$','fontsize',14)
ylabel('$F(pN)$','fontsize',14)
legend('$F_{tot}$','$F_2$')
end

if(0)   % FORZA F1 F2 E Ftot IN FUNZIONE DEL TEMPO
    figure(2)
plot(tempo/N, F(:,3),'-k'), hold on
plot(tempo/N, F(:,1),'Color', '0.64,0.08,0.18'), hold on
plot(tempo/N, F(:,2),'Color','0.00,0.45,0.74'), hold on
line([tempo(1)/N tempo(end)/N], [mean(Ftot) mean(Ftot)], 'Color', 'r')
xlabel('$t(s)$','fontsize',14)
ylabel('$F(pN)$','fontsize',14)
%xlim([0 0.32])
legend('$F_{tot}$','$F_1$','$F_2$', 'Fontsize', 10)
end

if(1)   % CONCENTRAZIONI IN FUNZIONE DEL TEMPO
    figure(4)
plot(tempo/N, nD/N, '-','Color', '#EDB120'),hold on
plot(tempo/N, nA1/N, '-','Color', '#A2142F'),hold on
plot(tempo/N, nA2/N, '-','Color','#0072BD'),hold on
plot(tempo/N, nD/N+nA1/N+nA2/N, '-','Color', 'black'   ,'Linewidth', 1),hold on
xlabel('$t\ (s)$','fontsize',14)
ylabel('$\displaystyle{\frac{n_D(t)}{N},\frac{n_1(t)}{N},\frac{n_2(t)}{N}}$','fontsize',14)
%xlim([0 0.35])
ylim([0 1.1])
legend('$D$','$A_1$','$A_2$')
end

if(0)   % CONCENTRAZIONI IN FUNZIONE DEL TEMPO (SUBPLOT)
    figure(4)
    subplot(3,1,1)
plot(tempo/N,nD/N, '-','Color', '#EDB120'), hold on        % 
%line([tempo(1)/N tempo(end)/N],[0.53 0.53], 'linestyle', '--'),hold on
%line([tempo(1)/N tempo(end)/N],[mean(nD/N) mean(nD/N)])
xlabel('$t\ s(s)$','fontsize',14)
ylabel('$\frac{n_D}{N}$','fontsize',14)
%xlim([0 0.35])
    subplot(3,1,2)
plot(tempo/N,nA1/N,'-', 'Color', '#A2142F'),hold on     %'Color','#77AC30'
%line([tempo(1)/N tempo(end)/N],[mean(nA1/N) mean(nA1/N)])
xlabel('$t\ (s)$','fontsize',14)
ylabel('$\frac{n_1}{N}$','fontsize',14)
%xlim([0 0.35])
    subplot(3,1,3)
plot(tempo/N,nA2/N,'-', 'Color', '#0072BD'),hold on     %'Color','#A2142F'
%line([tempo(1)/N tempo(end)/N],[mean(nA2/N) mean(nA2/N)])
xlabel('$t\ (s)$', 'fontsize',14)
ylabel('$\frac{n_2}{N}$','fontsize',14)
%xlim([0 0.35])

%    figure(6)
%hist(nA1/N), hold on
 %   figure(7)
%hist(nA2/N), hold on

end

if(0)   % ATPasi TOTALE IN FUNZIONE DEL TEMPO
    figure(5)
plot(tempo/N, atpTot/N,'o'), hold on
xlabel('$t(s)$', 'fontsize',14)
ylabel('$\phi (N^{-1} s^{-1})$','fontsize',14)
end

if(0)   % FLUTTUAZIONI Xi, Eta, XiEta IN FUNZIONE DEL TEMPO
  figure(1)
  plot(tempo/N, xi, '-', 'Color', '#EDB120'), hold on
  xlabel('$t$','fontsize',14)
  ylabel('$\xi_1$','fontsize',14)
    figure(2)
  plot(tempo/N, eta, '-', 'Color', '#EDB120'),hold on
  xlabel('$t$','fontsize',14)
  ylabel('$\xi_2$','fontsize',14)
    figure(3)
  plot(tempo/N, xieta, '-', 'Color', '#EDB120'),hold on
  xlabel('$t$','fontsize',14)
  ylabel('$\xi_1 \xi_2$','fontsize',14)
end

if (0)   % DISTRIBUZIONE DI PROBABILITA DELLA FLUTTUAZIONI
        figure(10)
    surf(xi, eta, P, 'LineWidth', 0.3), hold on
    xlabel('$\xi$', 'fontsize', 14)
    ylabel('$\eta$','fontsize', 14)
    zlabel('$\Pi(\xi, \eta)$', 'fontsize', 14)
        figure(11)
    plot(xi, Pxi, 'o-'), hold on
    xlabel('$\xi$', 'fontsize', 14)
    ylabel('$\Pi(\xi)$','fontsize', 14)
    [yxi xxi]=hist((xi-mean(xi)),5);
    yxi= yxi/sum(yxi)/(xxi(2)-xxi(1));
    plot(xxi, yxi), hold on
        figure(12)
    plot(eta, Peta, 'o-'), hold on
    xlabel('$\eta$', 'fontsize', 14)
    ylabel('$\Pi(\eta)$','fontsize', 14)
    [yeta xeta]=hist((eta-mean(eta)),5);
    yeta= yeta/sum(yeta)/(xeta(2)-xeta(1));
    plot(xeta, yeta), hold on
    
end  
%end

if(0)   % MOMENTI SECONDI DELLE FLUTTUAZIONI IN FUNCIONE DI K2
     %   %#EDB120 giallo     #64d413  verde #0072bd blu    #a2142f rosso
   figure(4)
  %plot(k2,Xi1_2, 'o', 'Color', 'r'), hold on
  %errorbar(k2,Xi11,errXi11, 'o', 'Color', 'm'), hold on  
  errorbar(storeZeta(:,1),storeZeta(:,2), storeZeta(:,3), 'o', 'Color', '#64d413', 'MarkerFaceColor', '#EDB120'), hold on
  xlabel('$k_2$','fontsize',14)
  ylabel('$<\xi^2>$','fontsize',14)
    figure(5)
  %plot(k2, Xi2_2, 'o', 'Color', 'r'),hold on
  %errorbar(k2,Xi22,errXi22, 'o', 'Color', 'm'), hold on
  errorbar(storeZeta(:,1),storeZeta(:,4), storeZeta(:,5), 'o', 'Color', '#64d413', 'MarkerFaceColor', '#EDB120'), hold on
  xlabel('$k_2$','fontsize',14)
  ylabel('$<\eta^2>$','fontsize',14)
    figure(6)
  %plot(k2, Xi1Xi2, 'o', 'Color', 'r'),hold on
  %errorbar(k2,Xi1Xi2,errXi12, 'o', 'Color', 'm'), hold on
  errorbar(storeZeta(:,1),storeZeta(:,6), storeZeta(:,7), 'o', 'Color', '#64d413', 'MarkerFaceColor', '#EDB120'), hold on
  xlabel('$k_2$','fontsize',14)
  ylabel('$<\xi \eta>$','fontsize',14)
end

if(1)
     figure(5)
plot(Xf,Yf/sum(Yf)/(Xf(2)-Xf(1)),'-k'), hold on
%bar(Xf,Yf/sum(Yf)/(Xf(2)-Xf(1)),'FaceColor','0.00,0.45,0.74','EdgeColor','k','LineWidth',0.5,'BarWidth',1), hold on
xlabel('$F(pN)$','fontsize',14)
ylabel('$P(F_{tot})$','fontsize',14)
    figure(6)
plot(Xf1,Yf1/sum(Yf1)/(Xf1(2)-Xf1(1))), hold on
%bar(Xf1,Yf1/sum(Yf1)/(Xf1(2)-Xf1(1)),'FaceColor','0.00,0.45,0.74','EdgeColor','k','LineWidth',0.5,'BarWidth',1), hold on
xlabel('$F(pN)$','fontsize',14)
ylabel('$P(F_1)$','fontsize',14)

    figure(7)
plot(Xf2,Yf2/sum(Yf2)/(Xf2(2)-Xf2(1))), hold on
bar(Xf2,Yf2/sum(Yf2)/(Xf2(2)-Xf2(1)),'FaceColor','0.00,0.45,0.74','EdgeColor','k','LineWidth',0.5,'BarWidth',1), hold on
xlabel('$F(pN)$','fontsize',14)
ylabel('$P(F_2)$','fontsize',14)

    figure(8)
plot(Xf, Yf, '-k'), hold on
plot(Xf1, Yf1, '-','Color','#77AC30'), hold on
plot(Xf2, Yf2, '-','Color','#A2142F'), hold on
%bar(Xf,Yf/sum(Yf)/(Xf(2)-Xf(1)),'FaceColor','0.00,0.45,0.74','EdgeColor','k','LineWidth',0.5,'BarWidth',1), hold on
xlabel('$F(pN)$','fontsize',14)
ylabel('$P(F_{tot}),\ P(F_1),\ P(F_2)$','fontsize',14)
xlim([-20 45])
legend('$P(F_{tot})$','$P(F_1)$','$P(F_2)$')

end

    figure(13)
plot(xConc, PD, ':k'),hold on
plot(xConc, P1, '-d', 'Color', '#f58453', 'MarkerEdgeColor', '#f58453', 'MarkerFaceColor', '#f58453'),hold on
plot(xConc, P2, '-o', 'Color', '#4dbfed', 'MarkerEdgeColor', '#4dbfed', 'MarkerFaceColor', '#4dbfed'),hold on
xlabel('$n_i/N$','fontsize',16)
ylabel('$PDF$','fontsize',14)
legend('$P(n_D/N)$','$P(n_1/N)$','$P(n_2/N)$')


if(saveResults)
    currDate = string(datetime, 'yyyy-MM-dd_HH-mm-ss');
    mkdir('datiGillespie',currDate)       
    outputName1= sprintf('datiGillespie/%s/datiGill_Popul_%s.dat',currDate,parSet);
    outputName2= sprintf('datiGillespie/%s/datiGill_Force_%s.dat',currDate,parSet);
 %   outputName3= sprintf('datiGillespie/%s/datiGill_DevF2Ftot_%s.dat',currDate,parSet);
    data1= Popul;
    data2= Force;
%    save(outputName3, '-ascii', 'data3');
    save(outputName1, '-ascii', 'data1');
    save(outputName2, '-ascii', 'data2');
    outputName4= sprintf('datiGillespie/%s/datiIstoGillPopA1_%s.dat',currDate,parSet);
    outputName5= sprintf('datiGillespie/%s/datiIstoGillPopA2_%s.dat',currDate,parSet);
    outputName6= sprintf('datiGillespie/%s/datiIstoGillPopD_%s.dat',currDate,parSet);  
    save(outputName4, '-ascii', 'P1');
    save(outputName5, '-ascii', 'P2');
    save(outputName6, '-ascii', 'PD');
    outputName7= sprintf('datiGillespie/%s/datiIstoGillForce_%s.dat',currDate,parSet);
    save(outputName7, '-ascii', 'istoF')
end
