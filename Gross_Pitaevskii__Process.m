%% Split-Operator Method
%% Ref: M. D. Feit, J. A. Fleck, Jr., and A. Steiger, 
%%      "Solution of the Gross-Pintaevskii Equation by a Spectral Method",
%%      Journal of Computational Physics 47, 412-433 (1982).
%% Ref: D. J. Griffiths, 
%%      "Introduction to Quantum Mechanics",
%%      ISBN 0-13-124405-1
%%  
%% See also: Gross_Pitaevskii__Process.m
%function so_teslavalve
%function y=Gross_Pitaevskii__Process(File_index,C_idx,P_idx)

function y=Gross_Pitaevskii__Process(dt,Adr,NPt)
tic
% v--- Change number in this section
%clear all; close all; clc;
%File_index = 103;  % Select file number as listed in outidxlist.m
%op = outidxlist(File_index);
%C_idx = 1;  % Select center number. 1 for left and 2 for right.
%P_idx = 2;  % Select the Mask number. 1 for left and 2 for right.
A_flag = 0; % Turn on/off gaussian absorption centers (sinks)
T_flag = 0; % Turn on/off contour of Tesla valve
Ndisp = 1;  %% Step interval to display
% NPt = 10000;%100000;%  %% Number of time step

op.center =[256,180;256,256;256,332];%[256,180; 500,500; 256,332];%[64,45; 64,64; 64,83];%[128,90; 128,128; 128,166]; %[256,180;256,256;256,332];%
op.dir = [1,0; -1,0];
op.M = 512;
    
%%
h = 6.6260693*10^-34; % [J s]
hbar = h/(2*pi); %[J s]
mRb = (87*(10^-3))/(6.0221415*(10^23)); % [kg]
Kb = 1.38*(10^-23); % [J K-1]
r0 = 5.2917721067e-11; %[m] Bohr radius
SLeng=100*r0;
g=0*4*pi*(hbar^2)*SLeng/mRb; %[coupling constant]

%% from /split_operator/bessel_scale.m
lambda_Rb=1e-6;
T0=((h/lambda_Rb)^2)/(2*mRb*Kb);
temp = 1e6*T0;%2e-6;%5e-7;%2e-6;%1e-5;% %% [K] Potential height
rscale = 1e-6;
rscaleinv = 1e6; %% rscaleinv = 1/rscale;
%dt = 1e-7;%1e-4;%1e-6;% 1e-5;% 1e-4;%  %% [s] Time step
%
%%
a = 10;  %% [rscale m] Length (Simulation dimension)
M = op.M;%512;%256;%128;
hM = M/2;
x = (rscale)*linspace(-a/2,a/2,M); x = x';
dx = x(2)-x(1);  % dx = a/M;
[xx,yy] = meshgrid(x,x);
% [theta,rho] = cart2pol(xx,yy);
%%
kx = linspace(-M/2,(M-2)/2,M); kx = kx';
% dkx = 1;
[kxx,kyy] = meshgrid(kx,kx);

%%-------------------------------------------------------------------------
%%
%%% Potential
%{
filename = op.name;
aV = imread(['D:\Work\Tesla_Valve\BECs\Analyse_Eigenenergy_Paper\Sub01\',filename]);
if File_index>=101
    %aV=im_resize(aV,.25); %Resize to match space
    aV=imresize(aV,0.25);
end

bV = double(aV(:,:,1));
%

V0 = ones(M,M);
xB = floor(size(bV,1)/2)-1;
yB = floor(size(bV,2)/2)-1;
V0(hM-xB:hM+1+xB+mod(size(bV,1),2),hM-yB:hM+1+yB+mod(size(bV,2),2)) = bV/255; %set up picture on a large plane
B = exp(-(((xx)/2e-6).^2 + ((yy)/2e-6).^2));
V0 = abs(ifft2c(fft2c(V0).*B));
figure(1);set(gcf,'position',[10 100 400 400]);
surf(x,x,V0);shading(gca,'interp');
% figure(2);set(gcf,'position',[10 500 400 400]);
% imagesc(bV);axis equal;axis tight;
%%
%}
V=[];
w=sqrt(2*pi/6.626*1e9);
V =(0.5)*(mRb)*w^2*(xx.^2+yy.^2);%0.5e-28*V0;%+Kb*temp*V0;
%V=ones(M,M)*5e-30;
%V(xx>-2e-6&xx<2e-6&yy>-2e-6&yy<2e-6)=0;
figure(2);set(gcf,'position',[400 100 500 400],'DefaultAxesFontSize',16);
surf(xx,yy,V); shading interp;
xlabel('x(m)'); ylabel('y(m)'); zlabel('Potential (J)');
%imagesc(V);axis equal;axis tight;axis xy;

%%-------------------------------------------------------------------------
%% Initial Gaussian wave packet (source)

Dy = (op.center(2,1)-hM)*dx;
Dx = (op.center(2,2)-hM)*dx;
xc = (rscale)*linspace(-a/2,a/2,M) - Dx; xc = xc';
yc = (rscale)*linspace(-a/2,a/2,M) - Dy; yc = yc';
[xxc,yyc] = meshgrid(xc,yc);
[theta,rhc] = cart2pol(xxc,yyc);
%Phi0 = exp( -1*0.5*(xxc.^2/4e-7) -1*0.5*(yyc.^2/4e-7)).*exp(1i*rhc);
%Phi0=2*exp(-1*rhc.^2/1e-12)+3*exp(-1*((xxc-2.5e-6).^2+(yyc).^2)/1e-12)+1*exp(-1*((xxc+2.5e-6).^2+(yyc).^2)/1e-12)+3*exp(-1*((xxc).^2+(yyc-2.5e-6).^2)/1e-12)+1*exp(-1*((xxc).^2+(yyc+2.5e-6).^2)/1e-12);
Phi0=exp(-1*0.5*rhc.^2/4e-10).*exp(1i*theta);
%Phi0=1./(1+rhc).*theta;

% + left->right up->down  - right->left down->up
%{
%%% Create Mask to calculate probability %%%
RPy = (op.center(3,1)-hM)*dx;
RPx = (op.center(3,2)-hM)*dx;
Rpxc = (rscale)*linspace(-a/2,a/2,M) - RPx; Rpxc = Rpxc';
Rpyc = (rscale)*linspace(-a/2,a/2,M) - RPy; Rpyc = Rpyc';
[xxpc,yypc] = meshgrid(Rpxc,Rpyc);
[~,Rprhc] = cart2pol(xxpc,yypc);
RPmask = Rprhc<0.5e-6; %positions where satisfy inequality are true, then return 1.
%%%%
%}
%%% Create Mask to calculate probability before output %%%
CPy = (op.center(2,1)-hM)*dx;
CPx = (op.center(2,2)-hM)*dx;
Cpxc = (rscale)*linspace(-a/2,a/2,M) - CPx; Cpxc = Cpxc';
Cpyc = (rscale)*linspace(-a/2,a/2,M) - CPy; Cpyc = Cpyc';
[xxCpc,yyCpc] = meshgrid(Cpxc,Cpyc);
[~,Cprhc] = cart2pol(xxCpc,yyCpc);
CPmask = Cprhc<Adr; %positions where satisfy inequality are true, then return 1.
%%%%'
%
%{
%%% Create Mask to absorb wavepacket %%%
Aby = (op.center(2,1)-hM+30)*dx;
Abx = (op.center(2,2)-hM)*dx;
Abxc = (rscale)*linspace(-a/2,a/2,M) - Abx; Abxc = Abxc';
Abyc = (rscale)*linspace(-a/2,a/2,M) - Aby; Abyc = Abyc';
[xxAbc,yyAbc] = meshgrid(Abxc,Abyc);
[~,Abrhc] = cart2pol(xxAbc,yyAbc);
Mask_absorb=zeros(M,M)+1;
Mask_absorb(Abrhc<0.5e-6)=0; %positions where satisfy inequality are true, then return 1.
Pmask_absorb=Abrhc<0.5e-6;
%}
%{
%%% Create Mask to initial wavepacket %%%
LPy = (op.center(2,1)-hM)*dx;
LPx = (op.center(2,2)-hM)*dx;
Lpxc = (rscale)*linspace(-a/2,a/2,M) - LPx; Lpxc = Lpxc';
Lpyc = (rscale)*linspace(-a/2,a/2,M) - LPy; Lpyc = Lpyc';
[xxLpc,yyLpc] = meshgrid(Lpxc,Lpyc);
[~,Lprhc] = cart2pol(xxLpc,yyLpc);
LPmask=Lprhc<5e-6;%0.5e-6;
%%%%
%}

Phi0 = Phi0/sqrt(sum(sum(abs(Phi0).^2)));
Phi0c = conj(Phi0); %% real(Phi0)- i*imag(Phi0);
%
figure(3);set(gcf,'position',[800 100 550 550]);
%subplot(2,1,1);
%contour(xx,yy,V,14);
xlabel('X(m)');ylabel('Y(m)');
%subplot(2,1,2);
contour(xx,yy,V,14);
hold on;
%contour(xx,yy,RPmask,14);
%contour(xx,yy,abs(Phi0),14);
%contour(xx,yy,Pmask_absorb,14);
contour(xx,yy,CPmask,14);
%contour(xx,yy,LPmask,14);
hold off;
pause(1);

figure(10);set(gcf,'position',[800 100 550 550]);
%subplot(2,1,1);
%contour(xx,yy,V,14);
%xlabel('X(m)');ylabel('Y(m)');
%subplot(2,1,2);
surf(xx,yy,V);

rectangle('Position',[0,0,2e-6,2e-6],'Curvature',[1 1])
shading interp;
%hold on;
%max(max(RPmask))
%imagesc(xx,yy,RPmask*0.5e-8);
%contour(xx,yy,abs(Phi0),14);
%contour(xx,yy,Pmask_absorb,14);
%imagesc(xx,yy,CPmask);
%imagesc(xx,yy,LPmask);
hold off;


figure(5);
set(gcf,'DefaultAxesFontSize',16);
TPhi0=Phi0./sum(sum(Phi0.*conj(Phi0)));
surf(xx,yy,real(TPhi0));
shading interp;
xlabel('x(m)'); ylabel('y(m)'); zlabel('Wavepacket Amplitude');
%ztickformat('%,.2f')

% Absorption around calculation boundary
Ax = exp(-((xx-x(1))/4e-7).^2)+exp(-((xx-x(end))/4e-7).^2);
Ay = exp(-((yy-x(1))/4e-7).^2)+exp(-((yy-x(end))/4e-7).^2);
A = Ax + Ay; A(A>1) = 1;
%
figure(4);set(gcf,'position',[10 50 400 400]);
imagesc(A);axis equal;axis tight;
%V = V -1i*Kb*temp*A*0;

Amask = zeros(M,M,3);

%% Sinks
if A_flag,
    Ac = op.center;  Ac(C_idx,:) = []; % remaining Absorption Centers
    figure(3); hold on;
    for nC = 1:size(Ac,1)
        Dy = (Ac(nC,1)-hM)*dx;
        Dx = (Ac(nC,2)-hM)*dx;
        xc = (rscale)*linspace(-a/2,a/2,M) - Dx; xc = xc';
        yc = (rscale)*linspace(-a/2,a/2,M) - Dy; yc = yc';
        [xxc,yyc] = meshgrid(xc,yc);
        [~,rhc] = cart2pol(xxc,yyc);
        Vsink = exp(-1*(rhc/2e-7).^2);
        Amask(:,:,nC) = rhc<0.6e-6;
        contour(xx,yy,Amask(:,:,nC),14);
        contour(xx,yy,abs(Vsink));
        V = V -1i*Kb*temp*Vsink;
    end
    figure(3); hold off;
end

%%-------------------------------------------------------------------------
%%
GK1 = (exp(-(1i*hbar*dt/(4*mRb))*((2*pi*rscaleinv/a)^2)*(kxx.^2 + kyy.^2))); 
% ^-- dt/2 kinetic energy propagator
GK2 = (exp(-(1i*hbar*dt/(2*mRb))*((2*pi*rscaleinv/a)^2)*(kxx.^2 + kyy.^2))); 
% ^-- dt kinetic energy propagator
GV0 = exp(-1i*(V+g*(abs(Phi0)).^2)*dt/hbar); %% Potential: spatial interaction
%%
iPhi = fft2c(Phi0);%*dx/sqrt(2*pi*hbar); %momentum space
figure(10);
surf(xx,yy,abs(iPhi)); shading interp; view(0,90);
Phi = ifft2c(iPhi.*GK1);%*sqrt(2*pi*hbar)/dx; %position space
Phi = GV0.*Phi;
%
Pt = zeros(1,NPt);
LPt = zeros(1,NPt);
ProbAll = zeros(1,NPt);
ProbMask = zeros(size(op.center,1),NPt);
ProbAbsorb =zeros(1,NPt);
ProbSMask=zeros(1,NPt);
%
% En = 2.276*1e-30;
%
%%
T = dt*NPt;
t = linspace(0,T,NPt);
Pt=zeros(1,NPt);
% wt = (1-cos(2*pi*t/length(t)));
% wt = (1-cos(2*pi*t/T));
% uns = 0;

figure(100);set(gcf,'position',[200 250 400 400]);
h100 = imagesc(x,x,abs(Phi0));hold on;
axis equal;axis tight;axis xy;
if T_flag,
    hV = contour(xx,yy,0.1*max(max(abs(Phi0)))*abs(V0),14);
end
title(['\Phi_x  t_{',num2str(1),'} = ',num2str(t(1))]);
xlabel('x (m)');ylabel('y (m)');

%Bk_Phi=zeros([size(Phi0),NPt]);
figure(100);set(gcf,'position',[900 200 400 400]);

%LPhi0=conj(LPmask.*Phi0);
%RPhi0=conj(RPmask.*Phi0);
CPhi0=conj(CPmask.*Phi0);
%need Shanks transformation
for nrn = 1:NPt
    iPhi = fft2c(Phi); %Momentum space
    Phi = ifft2c(iPhi.*GK1); %position space
    GV = exp(-1i*(V+g*(abs(Phi)).^2)*dt/hbar);
    PPhi = abs(Phi);
    %     Pt(nrn) = trapz(x,trapz(x,Phi0c.*Phi,2),1);
    %ProbAbsorb(1,nrn)=sum(sum((Pmask_absorb.*PPhi).^2)); %Prob that is absorped
    %Phi=Phi.*Mask_absorb; %Absorb wavepacket in specular area
    ProbSMask(1,nrn)=sum(sum((CPmask.*PPhi).^2)); %Prob that is still exist in tesla structure
    %Bk_Phi(:,:,nrn)=PPhi;
    %Normalization%
    
    %PPhi = abs(Phi)./sqrt(sum(sum(PPhi.^2)));
    %
    ProbAll(nrn) = sum(sum(PPhi.^2)); %Prob of all space
    %ProbMask(1,nrn) = sum(sum((Pmask.*PPhi).^2)); %Prob of output and before absorbtion area.
    
    %%Correlation of Left Center and Right area%%

   % LPt(nrn)=trapz(x,trapz(x,LPhi0.*LPmask*Phi));
    CPt(nrn)=trapz(x,trapz(x,CPhi0.*CPmask*Phi));
    %RPt(nrn)=trapz(x,trapz(x,RPhi0.*RPmask*Phi));
%Correlation Function is CPt
    
    if A_flag
        for nC = 1:size(Ac,1)
            ProbMask(nC+1,nrn) = sum(sum((Amask(:,:,nC).*PPhi).^2));
        end
    end
    %     unl = Phi*wt(nrn)*exp(1i*En*nrn*dt/hbar);
    %     if nrn > 1
    %         una = (unp + unl)*dt; %% Trapezoidal area
    %         uns = uns + una;
    %     end
    %     unp = unl;
   %{
    if mod(nrn,Ndisp) == 0
        figure(100);set(h100,'CData',PPhi);
        title(['\Phi_x  t_{',num2str(nrn),'} = ',num2str(t(nrn))]);
        contour(xx,yy,abs(Phi)); drawnow;%pause(0.05);
    end
    %}
    iPhi = fft2c(Phi);% Momentum space
    Phi = ifft2c(iPhi.*GK1); %Position space
    Phi = GV.*Phi;
    
end
%%
%Additional part
%Shanks transformation will be applied here
Shanks1=[];
%CPt(nrn) is square matrix
% disp('length of CPt is');disp(length(CPt));
% disp(CPt);%%the current problem is CPt has imaginary part
% for j=1:length(CPt)
%     Shanks=(CPt(j+2)*CPt(j)-(CPt(j+1))^2)/(CPt(j+2)-2*CPt(j+1)+CPt(j));%Shanks1
%     Shanks1=[Shanks1,Shanks];
% end
% N=length(CPt);%Shanks1);%CPt
% E=(1/N)*(1/dt)*linspace(-N/2,N/2-1,N)*h;
% CPe_sh = fftshift(fft(((1-cos(2*pi*t/T)).*CPt/T)));
% for j=1:length(CPe_sh)
%     Shanks=(CPe_sh(j+2)*CPe_sh(j)-(CPe_sh(j+1))^2)/(CPe_sh(j+2)-2*CPe_sh(j+1)+CPe_sh(j));%Shanks1
%     Shanks1=[Shanks1,Shanks];
% end
% figure(1004);%%
% hold on;
% plot(E,log(fliplr(abs(CPe_sh))),'k');%axis([0 3.5e-28 -30.5 -26.5]);
%end additional part
%%
% figure(101);set(gcf,'position',[400 150 800 500]);
% subplot(2,2,1);
% plot((1:NPt)*dt,ProbAll); axis([1*dt NPt*dt 0 1.5]);
% xlabel('t ({s})');ylabel('All Probability');
% %figure(102);set(gcf,'position',[600 200 400 400]);
% subplot(2,2,2);
% plot((1:NPt)*dt,ProbMask(1,:)); axis([1*dt NPt*dt 0 1]);
% xlabel('t {s})');ylabel('Output Probability'); axis tight;
% subplot(2,2,3);
% plot((1:NPt)*dt,ProbAbsorb(1,:)); axis([1*dt NPt*dt 0 1]);
% xlabel('t ({s})');ylabel('Absorb Probability'); axis tight;
% subplot(2,2,4);
% plot((1:NPt)*dt,ProbSMask(1,:)); axis([1*dt NPt*dt 0 1]);
% xlabel('t ({s})');ylabel('Input Probability'); axis tight;
% %suptitle('Probability vs time(s)');
% saveas(figure(101),['Probability_',num2str(File_index),'_',num2str(C_idx),'_',num2str(P_idx)],'tiffn');
% %}
N=length(CPt);
%E = (hbar/dt)*(linspace(-pi,pi,length(CPt)));
E=(1/N)*(1/dt)*linspace(-N/2,N/2-1,N)*h;
%%
%RPe = fftshift(fft(((1-cos(2*pi*t/T)).*RPt/T)));
CPe = fftshift(fft(((1-cos(2*pi*t/T)).*CPt/T)));
%LPe = fftshift(fft(((1-cos(2*pi*t/T)).*LPt/T)));
%Pe = fftshift(fft(Po));
%%
%Applied Shamks transformation
%for l=1:N
%S_A(l)=(A(l+2)*A(l)-(A(l+1))^2)/(A(l+2)-2*A(l+1)+A(l));

%%
%{
figure(102);subplot(2,1,1);plot(t,real(Pt));
title('Correlation Function ');xlabel('Time');
figure(102);subplot(2,1,2);plot(E,log(fliplr(abs(RPe))),'r');
title('Energy Spectrum');xlabel('Energy');ylabel('Power');
%}

%Pe = fftshift(fft(Po));
%%
%{
figure(103);subplot(2,1,1);plot(t,real(LPt));
title('Correlation Function ');xlabel('Time');
figure(103);subplot(2,1,2);plot(E,log(fliplr(abs(PeI))),'r');
title('Energy Spectrum');xlabel('Energy');ylabel('Power');
%}
figure(104);
%subplot(3,1,1);
hold on;
%plot(E,log(fliplr(abs(LPe))),'b');
%subplot(3,1,2);
 plot(E,log(fliplr(abs(CPe))),'k');%axis([0 3.5e-28 -30.5 -26.5]);
%subplot(3,1,3);
%plot(E,log(fliplr(abs(RPe))),'k');
%plot([0 0],[-45 -25],'--k');
%plot([100*Kb*T0 100*Kb*T0],[-45 -25],'--k');
%xlim([0 100*Kb*T0]);

K=[E;log(fliplr(abs(CPe)))];
title(['Energy Spectrum of Adr=2.8e-8M512upright']);xlabel('Energy');ylabel('Power');
legend('All space','Partial space at Center','Partial space at Right','Analytic Solution');
k=1;
num=50;
AnalyticSol=[];
for i=1:num
    for j=1:num
       % AnalyticSol(k)=((hbar*pi)^2 /(2*mRb))*((i/4e-6)^2+(j/0.4e-6)^2);
        AnalyticSol(k)=hbar*w*(0.5+0.5+i+j-2);
        k=k+1;
    end
end
stem(AnalyticSol,AnalyticSol*0-55,'Marker','*');
%stem(AnalyticSol,AnalyticSol*0-45,'Marker','*');
hold off;
%figure(104);saveas(gcf,['HarmonicOscillator_Adr_',num2str(Adr)],'fig');
%{
%Write average and SD of Probability
NonZeroProbMask=ProbMask;
NonZeroProbMask(NonZeroProbMask==0)=[]; % Calculate prob only non zero elements
NonZeroProbMask_average=mean(NonZeroProbMask);
NonZeroProbMask_SD=std(NonZeroProbMask);
fileID=fopen('Prob_Data_DoubleWell.txt','a+');% Non overwriting in text
fprintf(fileID,'%d %1d %1d %.5f %.5f\n',File_index,C_idx,P_idx,NonZeroProbMask_average,NonZeroProbMask_SD);
fclose(fileID);
%}
Processing_Time=toc;
saveas(gcf,'Adr0.028e-6M512upright.bmp')
%save(['HarmonicOscillator2D_Adr_',num2str(Adr)]);
%save(['CPe_Adr_',num2str(Adr),'.mat'],CPe,E);%close all;


%% Peaks finding %%
epC=log(fliplr(abs(CPe)));
pk=findpeaks(epC);
H=[];
%%find index of peaks
j=size(pk);
j=j(1,2);
%j=str2int(j);
for i=1:j  %%find the index in Y that are peaks
    G=find(epC==pk(i));
    H=[H,G];
end
I=[]; %%find the E that make peaks from the index of peaks of Y
for i=1:j
    k=H(i);
    n=E(k);
    I=[I,n];
end
[pkk,locx] = findpeaks(epC,E);
% display(locx);
%% find the min,max of AnalyticSol for restricting the interval of E by approximation
a=min(AnalyticSol);
b=max(AnalyticSol);
O1=find(a<=locx);
O2=find(locx<=b);
O=[];
if (numel(O2)>numel(O1)) | (numel(O2)==numel(O1)) 
    for i=1:numel(O2)
        O3=find(O1==O2(i));
        O=[O,O3];
    end
elseif (numel(O2)<numel(O1))
    for i=1:numel(O1)
        O3=find(O2==O1(i));
        O=[O,O3];
    end
end
A=find(I==a);
B=find(I==b);
Analytic_num=numel(AnalyticSol);
p=[];
display(O);
for i=min(O):max(O) %%restrict the interval of E
    c=I(i);
    p=[p,c];             %p is the E in the interval
end

% display('size p is');display(size(p));display(p); %display(I);
% findpeaks(I);
% display('size Analytic is');display(size(AnalyticSol));
% display(AnalyticSol);
AnSol=unique(AnalyticSol); %from AnalyticSol has repeated elements we have to cut it off
AnSole=AnSol(2:2:end); %for even level
AnSolo=AnSol(1:2:end); %ood
display('size unique Analytic is');display(size(AnSol));
%display(AnSol);
An_interv=[];
for i=1:numel(AnSol)-1  %find the interval of each level of Analytic
    in=AnSol(1,i+1)-AnSol(1,i);
    An_interv=[An_interv,in];
end
%display(An_interv);
%%For the E by approx has no unique peak in some levels%%
co=[];
App=[];
% if %???????????? Even Energy
	for i=1:numel(AnSole)	
                count1=find(I<(AnSole(i)+An_interv(i))); %find upper bound
                count2=find(I>(AnSole(i)-An_interv(i))); %find lower bound
                count=[min(count2):max(count1)];
                display(count);%display(count1);display(count2);
                display(I(min(count2)));display(I(max(count1)));display(i);
			
            if numel(count)==0
                
                
                if abs(AnSole(i)-I(min(count2)))>abs(AnSole(i)-I(max(count1)))
                    count=[max(count1)];
                elseif abs(AnSole(i)-I(min(count2)))<=abs(AnSole(i)-I(max(count1)))
                    count=[min(count2)];
                else
                end
            elseif count(1)==0
				count=count(2:end);
            else
			end
			coco=[];
            R=[];
			if numel(count)>0
				for k=min(count):max(count)
					coco1=I(k);
					coco=[coco,coco1];
				end
			R=max(coco);
            else
            end
        App=[App,R];
	end
% end	
    inter=An_interv;
if numel(App)==numel(AnSole)
    del=[(App-AnSole)./AnSole]*100;%abs
    %del=[];
    %for i=1:numel(App)
     %   del1=[abs(App(i)-AnSole(i))/(An_interv(i))]*100;
      %  del=[del,del1];
    %end
    figure(1000)
    elevel=2*find(AnSole);
    stem(elevel,del);
    title('Error between Analytic and self-correlation for Adr=0.028e-6M512upright');xlabel('level');ylabel('Error(%)');
    saveas(gcf,'ErrorAdr0.028e-6M512upright.bmp')
        % elseif numel(p)<Analytic_num %if missing a part of O
%     if p(1)>=AnalyticSol(2)
%         p=[I(A-1),p];
%     elseif (p(1)<AnalyticSol(2)) & (p(1)>((3*AnalyticSol(1)+AnalyticSol(2))/2))
%         p=[I(A-1),p];
% %     elseif O(1)<AnalyticSol(1)
% %         O=[I(A-1),O];
%     elseif p(numel(p))<((3*AnalyticSol(Analytic_num-1)+AnalyticSol(Analytic_num))/2)
%         p=[p,I(B+1)];
%     end
else
end
%%Saving data%%
Pe=[elevel.',del.']; %,AnalyticSol.',epC.',E.'
% eP=[AnalyticSol.',epC.',E.'];
save('test.mat','elevel','del')
dlmwrite('ErrorAdr0.028e-6M512upright.dat',Pe,'precision','%3.30e','delimiter',' ');
% dlmwrite('Adr0.2e-6.dat',eP,'precision','%3.30e','delimiter',' ');
% [pkk,locx] = findpeaks(epC,E);
% display(locx);

end


% display('H=');display(H);
% display('size of H');display(size(H));
% display('I=');display(I);
% display('size of I');display(size(I));



