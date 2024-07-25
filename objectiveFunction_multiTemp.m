function errorFunction = objectiveFunction_multiTemp(X,Data,DataT,DataDev,DataDevT,setPDF,DeltaT,Q1,Q_1,Q2,Q_2,Q3,N)

% Data to be fitted
tt= Data(:,1);
xx= Data(:,2);
ttT= DataT(:,1);
xxT= DataT(:,2);

ttDev= DataDev(:,1);
xxDev= DataDev(:,2);
ttDevT= DataDevT(:,1);
xxDevT= DataDevT(:,2);

% Fitting Parameters
f0= abs(X(1));
k1= abs(X(2));
k_1= abs(X(3));
k2= abs(X(4));
k_2= abs(X(5));
k3= abs(X(6));

k1T= k1*(Q1^(DeltaT/10));
k_1T= k_1*(Q_1^(DeltaT/10));
k2T= k2*(Q2^(DeltaT/10));
k_2T= k_2*(Q_2^(DeltaT/10));
k3T= k3*(Q3^(DeltaT/10));

[xm, xmDev]= pdfFitForce_2motors(k1,k_1,k2,k_2,k3,f0,N,Data,DataDev,setPDF);
[xmT, xmDevT]= pdfFitForce_2motors(k1T,k_1T,k2T,k_2T,k3T,f0,N,DataT,DataDevT,setPDF);

err1=(1/length(tt))*sum((xx-xm).^2);
err2= (1/length(ttT))*sum((xxT-xmT).^2);
errD= (1/length(ttDev))*sum((xxDev-xmDev).^2);
errDT= (1/length(ttDevT))*sum((xxDevT-xmDevT).^2);

errorFunction= err1+ err2+ (1/3)*(10^8)*errD+ (1/3)*(10^8)*errDT; %ME + Dev:MF 
%[err1 err2 (1/3)*(10^8)*errD (1/3)*(10^8)*errDT]

%errorFunction= err1+ err2+ 2*(10^-5)*errD+ 2*(10^-5)*errDT; %noise
%errorFunction= err1+ (10^-2)*err2+ 2*(10^-5)*errD+ 9*(10^-5)*errDT; %VK noise 34
%errorFunction= err1+ (10^-1)*err2+ (10^-6)*errD+ 2*(10^-6)*errDT; %VK noise

%errorFunction= err1+ err2+ (10^-5)*errD+ 0.5*(10^-5)*errDT; %ME + Dev noise 
%e=[err1 err2 (10^-5)*errD 0.5*(10^-5)*errDT];

%save errorfun e -ascii;

if(0)
    figure(10);
    hold off
    plot(tt,xx,'k');
    hold on
    plot(tt,xm,'r')
    figure(20);
    hold off
    plot(ttT,xxT,'k');
    hold on
    plot(ttT,xmT,'r')

    figure(30);
    hold off
    plot(ttDev,xxDev,'k');
    hold on
    plot(ttDev,xmDev,'r')

    figure(40);
    hold off
    plot(ttDevT,xxDevT,'k');
    hold on
    plot(ttDevT,xmDevT,'r')
end


