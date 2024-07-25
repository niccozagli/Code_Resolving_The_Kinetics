    % FIT DELLA SERIE TEMPORALE DELLA FORZA ISOMETRICA 2 MOTORI, 2 TEMPERATURE (ME/VK)

addpath('C:\Users\Vali\Desktop\MATLAB\PhD');
addpath(genpath('./'))
setFigure;

    % Selezionare il file dati e i parametri da fittare
expNum= '20240515';
parSet= 'SYNTH-A1A2-N20-f3-T25';
parSetT= 'SYNTH-A1A2-N20-f3-T15';
setPDF= 'ME'; 
%setPDF= 'VK';
setForceDev= 'MF';
%setForceDev= 'gillMedie'; 
saveResults= 1;
  
fileName= sprintf('../datiGillespie/%s/%s/datiIstoGillForce_%s.dat',expNum,parSet,parSet);
expName= [fileName(end-24:end-4) '-multiT'];
Data= load(fileName);
fMean= sum(Data(:,1).*Data(:,2))/sum(Data(:,2)) ;      %F media dell'ensemble

fileNameT= sprintf('../datiGillespie/%s/%s/datiIstoGillForce_%s.dat',expNum,parSetT,parSetT);
DataT= load(fileNameT);
fMeanT= sum(DataT(:,1).*DataT(:,2))/sum(DataT(:,2)) ;      %F media dell'ensemble
T0= str2double(parSet(end-1:end));
TT= str2double(parSetT(end-1:end));
DeltaT= TT-T0;

%taglia del sistema
N= 20;

switch setForceDev
    case 'gill'
        % %fileName= '../datiGillespie/SYNTH-A1A2_f03/datiIstoGillForce_SYNTH-A1A2-N50.dat';
        % %fileNameDev= '../datiGillespie/SYNTH-A1A2-f03/datiGill_Force_SYNTH-A1A2-N10000.dat';
        % Force= load(fileNameDev);
        % %DataDev= [Force(1:100000,1) Force(1:100000,4)];
        % DataDev= [Force(1:800,1) Force(1:800,4)];
        % %plot(DataDev(:,1), DataDev(:,2)),shg
    case 'gillMedie'
        fileNameDev= sprintf('../datiGillespie/%s/%s/datiGill_forceDev-av_%s_it800.dat',expNum,parSet,parSet);
        fileNameDevT= sprintf('../datiGillespie/%s/%s/datiGill_forceDev-av_%s_it800.dat',expNum,parSetT,parSetT);
       
        ForceDev= load(fileNameDev);
        %%Y= 0;%0.5*randn(length(ForceDev));
        %plot(DataDev(:,1), DataDev(:,2))
        DataDev= [ForceDev(:,1) ForceDev(:,2)];%0.05*rand(1)];
        %plot(DataDev(:,1), DataDev(:,2)),shg
        %DataDev= DataDev(1:300,:);
        DataDev= DataDev(1:4000,:);
        %plot(DataDev(:,1), DataDev(:,2)),shg

        ForceDevT= load(fileNameDevT);
        DataDevT= [ForceDevT(:,1) ForceDevT(:,2)];%0.05*rand(1)];
        DataDevT= DataDevT(1:4000,:);

        DataDev= [DataDev(:,1)/N DataDev(:,2)];  %quando si usano le medie va diviso il tempo per N
        DataDevT= [DataDevT(:,1)/N DataDevT(:,2)]; %quando si usano le medie va diviso il tempo per N
    case 'MF'
        %fileNameDev= '../datiMF_2motors/dataDev_MF_FtotDev-SYNTH-A1A2-N20-f3.dat';
        %fileNameDevT= '../datiMF_2motors/dataDev_MF_FtotDev-SYNTH-A1A2-N20-f3_T34.dat';
        
        fileNameDev= '../datiMF_2motors/dataDev_MF_SYNTH-A1A2-N20-f3-T25.dat';
        fileNameDevT= '../datiMF_2motors/dataDev_MF_SYNTH-A1A2-N20-f3-T15.dat';

        DataDev=load(fileNameDev);
        DataDevT=load(fileNameDevT);
end

    %fattori Q10 
Q1= 2;
Q_1= 2; 
Q2= 5.5;
Q_2= 1.5; 
Q3= 4;

    %FVAL max accettato
TOLL= 1*10^(-5);% ME+noise % 10^-5;%7*10^-7;%9*10^-7;% VK noise err1+ err2+ (10^-5)*errD+ 2*(10^-5)*errDT =2.9*10^-5;%VK MF err1+ err2+ (10^8)*errD+ (10^8)*errDT %3*10^-5; %e= err1+ err2+ (10^8)*errD+ (10^8)*errDT; %ME con Dev:MF;%4*10^-5;%3*10^-6;%11*10^-6; %ME

parTrue= [3 5 15 150 10 3];
%lb= [3 5 15 150 10 3];
%ub=[3 5 15 150 10 3];
lb= [0 0 0 0 0 0];
ub= [10 200 200 2000 500 200];

    %numero iterazioni
iteraz= 10;   
    %cond iniziali del fit
%startPar= [10 250 250 1000 250 250].*rand(iteraz,length(lb));
startPar= ub.*rand(iteraz,length(lb));
%startPar= parTrue;

V= [];
storeFit= [];  
storeFDev= [];
%ParBestStore= [];
sp= [];

for it= 1:iteraz
    [N it]    
    
    %P_subset= parTrue;
    %P_subset= lb+(ub-lb).*rand(1,length(lb));
    %P_subset= (ub./100).*rand(1,length(lb));
    P_subset= startPar(it,:);
  
    %OPTIONS = SIMPSASET('MAX_ITER_TOTAL',50000,'TOLFUN',1e-16, 'TOLX',1e-30);
    OPTIONS = SIMPSASET('TOLFUN', 1e-6, 'MAX_TIME',1000, 'MAX_FUN_EVALS',1000);
    tic
    [ParBest,FVAL,EXITFLAG,OUTPUT] = SIMPSA('objectiveFunction_multiTemp',P_subset,lb,ub,OPTIONS,Data,DataT,DataDev,DataDevT,setPDF,DeltaT,Q1,Q_1,Q2,Q_2,Q3,N);
    %[ParBest,FVAL,EXITFLAG,OUTPUT] = SIMPSA('objectiveFunction_multiTemp',P_subset,lb,ub,[],Data,DataT,DataDev,DataDevT,setPDF,DeltaT,Q1,Q_1,Q2,Q_2,Q3,N);
    if EXITFLAG ~=1
        sprintf('EXITFLAG=%d, it=%d', EXITFLAG,it)
    end
   toc
    f0= ParBest(1);
    k1= ParBest(2);
    k_1= ParBest(3);
    k2= ParBest(4);
    k_2= ParBest(5);
    k3= ParBest(6);
    G= (k_1*(k_2+k3)+k2*k3)/(k2+k_2+k3);
    r= k1/(k1+G);
    phi= (k1*k2*k3)/((k1+G)*(k2+k_2+k3));

    [PF, FDev]= pdfFitForce_2motors(k1,k_1,k2,k_2,k3,f0,N,Data,DataDev,setPDF);
    PF= PF';
   
    k1T= k1*(Q1^(DeltaT/10));
    k_1T= k_1*(Q_1^(DeltaT/10));
    k2T= k2*(Q2^(DeltaT/10));
    k_2T= k_2*(Q_2^(DeltaT/10));
    k3T= k3*(Q3^(DeltaT/10));

    GT= (k_1T*(k_2T+k3T)+k2T*k3T)/(k2T+k_2T+k3T);
    rT= k1T/(k1T+GT);
    phiT= (k1T*k2T*k3T)/((k1T+GT)*(k2T+k_2T+k3T));

    [PFT, FDevT]= pdfFitForce_2motors(k1T,k_1T,k2T,k_2T,k3T,f0,N,DataT,DataDevT,setPDF);
    PFT= PFT';
  
    if FVAL<TOLL
        V= [V
            N f0 k1 k_1 k2 k_2 k3 r rT phi phiT FVAL fMean fMeanT]
        sp=[sp
            startPar(it,:)];
        storeFit= [storeFit
            N PF PFT];
        storeFDev= [storeFDev
            N FDev' FDevT'];
    end
end


if(0)

    figure(1)
    bar(DataT(:,1),DataT(:,2),'Facealpha', 0.5, 'FaceColor','0.00,0.45,0.74','EdgeColor','k','LineWidth',0.5,'BarWidth',1),hold on
    plot(DataT(:,1),PFT,'Color','0.87,0.46,0.29','LineWidth',1),shg
    bar(Data(:,1),Data(:,2),'FaceColor','0.00,0.45,0.74','EdgeColor','k','LineWidth',0.5,'BarWidth',1),hold on
    plot(Data(:,1),PF,'Color','0.87,0.33,0.12','LineWidth',1),shg
    legend('$P(F,T_2)$', 'best fit ($T_2$)', '$P(F,T_1)$', 'best fit ($T_1$)')
    xlabel('$F$','fontsize',16)
    ylabel('$P(F)$','fontsize',16)

    figure(2)
    plot(DataDev(1:length(FDev),1),DataDev(1:length(FDev),2),'LineWidth',2),hold on
    plot(DataDev(1:length(FDev),1),FDev,'Color','0.85,0.33,0.10','LineWidth',1), shg
    plot(DataDevT(1:length(FDevT),1),DataDevT(1:length(FDevT),2),'LineWidth',2),hold on
    plot(DataDevT(1:length(FDevT),1),FDevT,'Color','0.85,0.33,0.10','LineWidth',1), shg
    legend('$F(t, T_1)$', 'best fit ($T_1$)', '$F(t, T_2)$', 'best fit ($T_2$)', 'location', 'southeast');
    ylabel('$F(t)$','fontsize',16)
    xlabel('$t$','fontsize',16)

    %         figure(1)
% bar(DataT(:,1),DataT(:,2),'Facealpha', 0.5, 'FaceColor','0.00,0.45,0.74','EdgeColor','k','LineWidth',0.5,'BarWidth',1),hold on
% plot(DataT(:,1),storeFit(:,32:end)','Color','0.87,0.46,0.29','LineWidth',1),hold on
% bar(Data(:,1),Data(:,2),'FaceColor','0.00,0.45,0.74','EdgeColor','k','LineWidth',0.5,'BarWidth',1),hold on
% plot(Data(:,1),storeFit(:,2:31)','Color','0.87,0.33,0.12','LineWidth',1),hold on
% legend('$P(F,T_2)$', 'best fit ($T_2$)', '$P(F,T_1)$', 'best fit ($T_1$)')
% xlabel('$F$','fontsize',16)
% ylabel('$P(F)$','fontsize',16)
% figure(2)
% plot(DataDev(:,1),DataDev(:,2),'LineWidth',2),hold on
% plot(DataDev(:,1),storeFDev(:,2:size(DataDev)+1),'Color','0.85,0.33,0.10','LineWidth',1), shg
% plot(DataDevT(:,1),DataDevT(:,2),'LineWidth',2),hold on
% plot(DataDevT(:,1),storeFDev(:,size(DataDev)+2:end),'Color','0.85,0.33,0.10','LineWidth',1), shg
% legend('$F(t, T_1)$', 'best fit ($T_1$)', '$F(t, T_2)$', 'best fit ($T_2$)', 'location', 'southeast');
% ylabel('$F(t)$','fontsize',16)
% xlabel('$t$','fontsize',16)

KT= V(:,3:7);

h3= figure(3);
hold on
plot(KT','.')
plot(parTrue(2:end)','or')
plot([1.1 2.1 3.1 4.1 5.1], startPar(:,2:end)', '.y');
errorbar(mean(KT), std(KT),'*k')
plot(lb(2:end),'xr')
plot(ub(2:end),'xr')
xlim([0 6])
xticks([1 2 3 4 5])
xticklabels({'$k_1$', '$k_{-1}$','$k_2$', '$k_{-2}$', '$k_3$'})

h4= figure(4);
violinplot(V(:,3:7)), hold on
plot(parTrue(2:end)','or')
end

%title('$T2=24^{\circ} \ C$')
%  %   case 2
%         KLT=V(:,11:15);
%         [p1LT p_1LT p2LT p_2LT p3LT]= readReactionConst('SYNTH-T10-N20-f3');
%         h4= figure(4);
%         plot(KLT','.'),hold on
%         plot([p1LT p_1LT p2LT p_2LT p3LT]','or'),hold on
%         errorbar(mean(KLT), std(KLT),'*k'),hold on
%         %plot(lb,'xr'),hold on
%         %plot(ub,'xr'),hold on
%         xlim([0 6])
%         xticks([1 2 3 4 5])
%         xticklabels({'$k_1$', '$k_{-1}$','$k_2$', '$k_{-2}$', '$k_3$'})
%         title('$LT$')
%         %title('$TLT=12^{\circ} \ C$')
%         %set(gca, 'YScale','log')
%         %legend('expected value', 'data 1', 'data 2', 'data 3', 'data 4', 'data 5','data 6', '$(mean, std)$', '$max$', '$min$', 'location', 'southeast')
% 
%     case 3
%         KT=V(:,3:7);
%         [p1 p_1 p2 p_2 p3]= readReactionConst('SYNTH-T24-N20-f3');
%         h3= figure(3);
%         plot(KT','.'),hold on
%         plot([p1 p_1 p2 p_2 p3]','or'),hold on
%         errorbar(mean(KT), std(KT),'*k'),hold on
%         %plot(lb,'xr'),hold on
%         %plot(ub,'xr'),hold on
%         xlim([0 6])
%         xticks([1 2 3 4 5])
%         xticklabels({'$k_1$', '$k_{-1}$','$k_2$', '$k_{-2}$', '$k_3$'})
%         title('$T2=24^{\circ} \ C$')
%         % KLT= V(:,11:15);
%         % KHT= V(:,16:20);
%         % [p1LT p_1LT p2LT p_2LT p3LT]= readReactionConst('SYNTH-T10-N20-f3');
%         % [p1HT p_1HT p2HT p_2HT p3HT]= readReactionConst('SYNTH-T32-N20-f3');
%         % h4= figure(4);
%         % plot(KLT','.'),hold on
%         % plot([p1LT p_1LT p2LT p_2LT p3LT]','or'),hold on
%         % errorbar(mean(KLT), std(KLT),'*k'),hold on
%         % %plot(lb,'xr'),hold on
%         % %plot(ub,'xr'),hold on
%         % xlim([0 6])
%         % xticks([1 2 3 4 5])
%         % xticklabels({'$k_1$', '$k_{-1}$','$k_2$', '$k_{-2}$', '$k_3$'})
%         % title('$LT$')
%         % %title('$TLT=12^{\circ} \ C$')
%         % %set(gca, 'YScale','log')
%         % %legend('expected value', 'data 1', 'data 2', 'data 3', 'data 4', 'data 5','data 6', '$(mean, std)$', '$max$', '$min$', 'location', 'southeast')
%         % h5= figure(5);
%         % plot(KHT','.'),hold on
%         % plot([p1HT p_1HT p2HT p_2HT p3HT]','or'),hold on
%         % errorbar(mean(KHT), std(KHT),'*k'),hold on
%         % %plot(lb,'xr'),hold on
%         % %plot(ub,'xr'),hold on
%         % xlim([0 6])
%         % xticks([1 2 3 4 5])
%         % xticklabels({'$k_1$', '$k_{-1}$','$k_2$', '$k_{-2}$', '$k_3$'})
%         % title('$HT$')
%         % %title('$T3=32^{\circ} \ C$')
%         % %set(gca, 'YScale','log')
%         % %legend('expected value', 'data 1', 'data 2', 'data 3', 'data 4', 'data 5','data 6', '$(mean, std)$', '$max$', '$min$', 'location', 'southeast')
% end


    % Si salvano i dati del fit e i valori dei parametri
if(saveResults)
    currDate = string(datetime, 'yyyy-MM-dd_HH-mm-ss');
    dirName= [currDate];
    mkdir('FiguresFit/2motors/fitCompleto', dirName)
    outputName= sprintf('FiguresFit/2motors/fitCompleto/%s/fitData_%s_%d-%d.dat',dirName,expName,T0,TT);
    outputNameDev= sprintf('FiguresFit/2motors/fitCompleto/%s/fitDataDev_%s_%d-%d.dat',dirName,expName,T0,TT);
    outputNamePar= sprintf('FiguresFit/2motors/fitCompleto/%s/fitParBest_%s_%d-%d.dat',dirName,expName,T0,TT);
    fitParBest= [V];
    save(outputNamePar, '-ascii', 'fitParBest');
    data3= [storeFit];
    data4= [storeFDev];
    save(outputName, '-ascii', 'data3');
    save(outputNameDev, '-ascii', 'data4');
    outputNameStart= sprintf('FiguresFit/2motors/fitCompleto/%s/startPar_%s_%d-%d.dat',dirName,expName,T0,TT);
    outputNameSP= sprintf('FiguresFit/2motors/fitCompleto/%s/sP_%s_%d-%d.dat',dirName,expName,T0,TT);
    save(outputNameStart, '-ascii', 'startPar');
    save(outputNameSP, '-ascii', 'sp');
end
