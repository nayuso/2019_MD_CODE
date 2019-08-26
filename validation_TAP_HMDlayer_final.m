%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VALIDATION OF THE ANALYTICAL EXPRESSIONS DUE TO A HMD
% IN A THREE-LAYERED REGION BY MEANS OF FEM SIMULATIONS
%
% Covered range: From near field to far field
%
%
% validation_HMDlayer.m
%
%
% Author: Dr. Natalia Ayuso
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
close all
clear all

global mu0 e0 c

computations=0;     %1 for calculus and 0 for loading stored computations
saving=0;           %1 for saving computations and 0 for avoiding saving

e0=8.85e-12;
mu0=4*pi*1e-7;
c=1/sqrt(e0*mu0);

% Emitter
m0=1;
he=-7;
f=[1e3 100e3 1e7];

% Media parameters
sigma1=10e-3;
er1=4;
e1=er1*e0;
sigma2=1e-3;
er2=4;
e2=er2*e0;
d=5;

% Receiver
hr_vector=[-2];
rx_hr=length(hr_vector);
rho=[0:1:100];
py=rho;
px=zeros(1,length(py));
rx_hr=length(hr_vector);
rx_rho=length(rho);


%MATLAB SIMULATIONS
dipo=1;         %dipole approximation
expresionx=0;   %Hx
expresiony=1;   %Hy
expresionz=2;   %Hz

aprox=0;
TOL=1E-15;

H_MAT1k_axis_xy=[];
H_MAT1k_axis_y=[];
H_MAT1k_axis_y2=[];

H_MAT100k_axis_xy=[];
H_MAT100k_axis_y=[];
H_MAT100k_axis_y2=[];

H_MAT10M_axis_xy=[];
H_MAT10M_axis_y=[];
H_MAT10M_axis_y2=[];


    
fprintf('Emitter placed in region 2 \n');
load( '../../MATLAB_FEM_DATA/HMDlayered_soil2_MATLAB_fine');
load( '../../MATLAB_FEM_DATA/HMDlayered_soil2_MATLAB_fine_LIM40_TOL14');
solHF=1;
 
fprintf('Field computation time: %3.0f min %2.0f s\n',minutes,segs);

if solHF==1
        fprintf('High Frecuency Solution\n')
        pause
end

Hx1kMAT_xy= abs(Hx_1k_MAT_axis_xy);     %Hx field along the xy axis
Hy1kMAT_y= abs(Hy_1k_MAT_axis_y);        %Hy field along the y axis
Hz1kMAT_y2= abs(Hz_1k_MAT_axis_y2);      %Hz field along the y axis but -0.5 m for he
Hx100kMAT_xy= abs(Hx_100k_MAT_axis_xy);
Hy100kMAT_y= abs(Hy_100k_MAT_axis_y);
Hz100kMAT_y2= abs(Hz_100k_MAT_axis_y2);
Hx10MMAT_xy= abs(Hx_10M_MAT_axis_xy);
Hy10MMAT_y= abs(Hy_10M_MAT_axis_y);
Hz10MMAT_y2= abs(Hz_10M_MAT_axis_y2);

style1k=['-+b -.+b--+b'];       %A line style for any receiver heigth and frequency
style100k=['-*g -.*g--*g'];     %A line style for any receiver heigth and frequency
style10M=['-or -.or--or'];      %A line style for any receiver heigth and frequency

fprintf('Solid line for hr=1, dotted for hr=-2 and dashed for hr=-7\n')

style1k=['-+b -.+b--+b']; %A line style for any receiver heigth and frequency
style100k=['-*g -.*g--*g']; %A line style for any receiver heigth and frequency
style10M=['-or -.or--or']; %A line style for any receiver heigth and frequency

fprintf('Solid line for hr=1, dotted for hr=-2 and dashed for hr=-7\n')

figure(1)
for i=1:rx_hr
    semilogy(R(i,:),Hx1kMAT_xy(i,:),style1k([(i-1)*4+1:i*4]),'LineWidth',2)
    hold on
    semilogy(R(i,:),Hx100kMAT_xy(i,:),style100k([(i-1)*4+1:i*4]),'LineWidth',2)
    semilogy(R(i,:),Hx10MMAT_xy(i,:),style10M([(i-1)*4+1:i*4]),'LineWidth',2)
end
set(gca,'FontSize',12,'Fontname','arial')
xlabel('R')
grid on
ylabel('|H_x|')
title('MATLAB SOLUTION along the xy axis')
legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])

figure(2)
for i=1:rx_hr
    semilogy(R(i,:),Hy1kMAT_y(i,:),style1k([(i-1)*4+1:i*4]),'LineWidth',2)
    hold on
    semilogy(R(i,:),Hy100kMAT_y(i,:),style100k([(i-1)*4+1:i*4]),'LineWidth',2)
    semilogy(R(i,:),Hy10MMAT_y(i,:),style10M([(i-1)*4+1:i*4]),'LineWidth',2)
end
set(gca,'FontSize',12,'Fontname','arial')
xlabel('R')
grid on
ylabel('|H_y|')
title('MATLAB SOLUTION along the y axis')
legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])

figure(3)
for i=1:rx_hr
    semilogy(Rz(i,:),Hz1kMAT_y2(i,:),style1k([(i-1)*4+1:i*4]),'LineWidth',2)
    hold on
    semilogy(Rz(i,:),Hz100kMAT_y2(i,:),style100k([(i-1)*4+1:i*4]),'LineWidth',2)
    semilogy(Rz(i,:),Hz10MMAT_y2(i,:),style10M([(i-1)*4+1:i*4]),'LineWidth',2)
end
set(gca,'FontSize',12,'Fontname','arial')
xlabel('R')
grid on
ylabel('|H_z|')
title('MATLAB SOLUTION along the y axis -0.5 elevation at the Tx')
legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])

%IMPORT COMSOL DATA. Don't forguet to write down the file extention
fprintf('Emitter placed in region 2 \n');
comsol_sol1kx=load('../../MATLAB_FEM_DATA/Hx_1k_med2_800.txt');
comsol_sol1ky=load('../../MATLAB_FEM_DATA/Hy_1k_med2_800.txt');
comsol_sol1kz=load('../../MATLAB_FEM_DATA/Hz_1k_med2_1600_ref2.txt');
comsol_sol100kx=load('../../MATLAB_FEM_DATA/Hx_100k_med2_final_ref.txt');
comsol_sol100ky=load('../../MATLAB_FEM_DATA/Hy_100k_med2_final_ref.txt');
comsol_sol100kz=load('../../MATLAB_FEM_DATA/Hz_100k_med2_final_ref.txt');      
comsol_sol10Mx=load('../../MATLAB_FEM_DATA/Hx_10M_med2_100.txt');
comsol_sol10My=load('../../MATLAB_FEM_DATA/Hy_10M_soil2_test.txt');
comsol_sol10Mz=load('../../MATLAB_FEM_DATA/Hz_10M_soil2_test.txt');
   

%Notes:  very importan. Due to simmetry the field is doubled. Consequently,
% field valueas have to be halved.
    
Hx_comsol=[comsol_sol1kx(:,4) comsol_sol100kx(:,4) comsol_sol10Mx(:,4)]/2;
Hy_comsol=[comsol_sol1ky(:,4) comsol_sol100ky(:,4) comsol_sol10My(:,4)]/2;
Hz_comsol=[comsol_sol1kz(:,4) comsol_sol100kz(:,4) comsol_sol10Mz(:,4)]/2;
  
rx_hr_comsol=3;   %Receiving heights
rx_rho_comsol=101; %Receiving points

style_comsol=['- -.--'];

figure(4)
%Hx
x_comsol=comsol_sol1kx(:,1);
y_comsol=comsol_sol1kx(:,2);
z_comsol=comsol_sol1kx(:,3);
for i=1:rx_hr_comsol
    R_comsol=sqrt(x_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+y_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(z_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol)-he).^2);
    semilogy(R_comsol,abs(Hx_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,:)),style_comsol([(i-1)*2+1:i*2]))
    hold on
end
set(gca,'FontSize',12,'Fontname','arial')
xlabel('R')
grid on
ylabel('|H_x|')
title('COMSOL SOLUTION along the xy axis')
legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])
    
figure(5)
%Hy
x_comsol=comsol_sol1ky(:,1);
y_comsol=comsol_sol1ky(:,2);
z_comsol=comsol_sol1ky(:,3);
for i=1:rx_hr_comsol
    R_comsol=sqrt(x_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+y_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(z_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol)-he).^2);
    semilogy(R_comsol,abs(Hy_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,:)),style_comsol([(i-1)*2+1:i*2]))
    hold on
end
set(gca,'FontSize',12,'Fontname','arial')
xlabel('R')
grid on
ylabel('|H_y|')
title('COMSOL SOLUTION along the y axis')
legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])
   
figure(6)
%Hz
x_comsol=comsol_sol1kz(:,1);
y_comsol=comsol_sol1kz(:,2);
z_comsol=comsol_sol1kz(:,3);
for i=1:rx_hr_comsol
    R_comsol=sqrt(x_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+y_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(z_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol)-he).^2);
    semilogy(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,:)),style_comsol([(i-1)*2+1:i*2]))
    hold on
end
set(gca,'FontSize',12,'Fontname','arial')
xlabel('R')
grid on
ylabel('|H_z|')
title('COMSOL SOLUTION along the y axis -0.5 elevation at the Tx')
legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])

figure(7)
%Hx
for i=1:rx_hr
    loglog(R(i,:),Hx1kMAT_xy(i,:),'+b')
    hold on
    loglog(R(i,:),Hx100kMAT_xy(i,:),'sg')
    loglog(R(i,:),Hx10MMAT_xy(i,:),'or')
end
x_comsol=comsol_sol1kx(:,1);
y_comsol=comsol_sol1kx(:,2);
z_comsol=comsol_sol1kx(:,3);
for i=1:rx_hr_comsol
    R_comsol=sqrt(x_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+y_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(z_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol)-he).^2);
    loglog(R_comsol,abs(Hx_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,:)),style_comsol([(i-1)*2+1:i*2]))
    hold on
end
set(gca,'FontSize',12,'Fontname','arial')
xlabel('R')
xlim([0 100])
grid on
ylabel('|H_x|')
title('MATLAB COMSOL SOLUTION')
legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])

figure(8)
%Hy
for i=1:rx_hr
    loglog(R(i,:),Hy1kMAT_y(i,:),'+b')
    hold on
    loglog(R(i,:),Hy100kMAT_y(i,:),'sg')
    loglog(R(i,:),Hy10MMAT_y(i,:),'or')
end
x_comsol=comsol_sol1ky(:,1);
y_comsol=comsol_sol1ky(:,2);
z_comsol=comsol_sol1ky(:,3);
for i=1:rx_hr_comsol
    R_comsol=sqrt(x_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+y_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(z_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol)-he).^2);
    loglog(R_comsol,abs(Hy_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,:)),style_comsol([(i-1)*2+1:i*2]))
    hold on
end
set(gca,'FontSize',12,'Fontname','arial')
xlabel('R')
xlim([0 100])
grid on
ylabel('|H_y|')
title('MATLAB COMSOL SOLUTION')
legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])

figure(9)
for i=1:rx_hr
    loglog(Rz(i,:),Hz1kMAT_y2(i,:),'+b')
    hold on
    loglog(Rz(i,:),Hz100kMAT_y2(i,:),'sg')
    loglog(Rz(i,:),Hz10MMAT_y2(i,:),'or')
end
x_comsol=comsol_sol1kz(:,1);
y_comsol=comsol_sol1kz(:,2);
z_comsol=comsol_sol1kz(:,3);
for i=1:rx_hr_comsol
    R_comsol=sqrt(x_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+y_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(z_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol)-he).^2);
    loglog(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,:)),style_comsol([(i-1)*2+1:i*2]))
    hold on
end
set(gca,'FontSize',12,'Fontname','arial')
xlabel('R')
xlim([0 100])
grid on
ylabel('|H_z|')
title('MATLAB COMSOL SOLUTION')
legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])

stepMAT=[2 3:2:101];
first=2;
rangeRcomsol=[first:length(R_comsol)]; %To avoid the first value that is NaN in numerical integration

%Tx=1 Rx=-2
 figure(10)
 %Hx
 i=2
 semilogy(R(i,stepMAT),Hx1kMAT_xy(i,stepMAT),'sb')
 hold on
    semilogy(R(i,stepMAT),Hx100kMAT_xy(i,stepMAT),'dk')
    semilogy(R(i,stepMAT),Hx10MMAT_xy(i,stepMAT),'or')
    x_comsol=comsol_sol1kx(:,1);
    y_comsol=comsol_sol1kx(:,2);
    z_comsol=comsol_sol1kx(:,3);
    R_comsol=sqrt(x_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+y_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(z_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol)-he).^2);
    semilogy(R_comsol(rangeRcomsol),abs(Hx_comsol((i-1)*rx_rho_comsol+first:i*rx_rho_comsol,1)),'-b','Linewidth',1.5)
    semilogy(R_comsol(rangeRcomsol),abs(Hx_comsol((i-1)*rx_rho_comsol+first:i*rx_rho_comsol,2)),':k','Linewidth',1.5)
    semilogy(R_comsol(rangeRcomsol),abs(Hx_comsol((i-1)*rx_rho_comsol+first:i*rx_rho_comsol,3)),'--r','Linewidth',1.5)

     set(gca,'FontSize',12,'Fontname','arial')
     xlabel('R (m)')
     xlim([0 101])
     ylim([1e-10 1])
     %grid on
     ylabel('|H_x| (A/m)')
     %legend(['f=1 kHz Analytic'],['f=100 kHz Analytic'],['f=10 MHz Analytic'], ['f=1kHz FEM'], ['f=100 kHz FEM'], ['f=10 MHz FEM'])
     legend(['f=1 kHz Numerical integration'],['f=100 kHz Numerical integration'],['f=10 MHz Numerical integration'], ['f=1kHz FEM computation'], ['f=100 kHz FEM computation'], ['f=10 MHz FEM computation'])

    figure(11)
    %Hy
    semilogy(R(i,stepMAT),Hy1kMAT_y(i,stepMAT),'sb')
    hold on
    semilogy(R(i,stepMAT),Hy100kMAT_y(i,stepMAT),'dk')
    semilogy(R(i,stepMAT),Hy10MMAT_y(i,stepMAT),'or')
    x_comsol=comsol_sol1ky(:,1);
    y_comsol=comsol_sol1ky(:,2);
    z_comsol=comsol_sol1ky(:,3);
    R_comsol=sqrt(x_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+y_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(z_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol)-he).^2);
    semilogy(R_comsol(rangeRcomsol),abs(Hy_comsol((i-1)*rx_rho_comsol+first:i*rx_rho_comsol,1)),'-b','Linewidth',1.5)
    semilogy(R_comsol(rangeRcomsol),abs(Hy_comsol((i-1)*rx_rho_comsol+first:i*rx_rho_comsol,2)),':k','Linewidth',1.5)
    semilogy(R_comsol(rangeRcomsol),abs(Hy_comsol((i-1)*rx_rho_comsol+first:i*rx_rho_comsol,3)),'--r','Linewidth',1.5)

     set(gca,'FontSize',12,'Fontname','arial')
     xlabel('R (m)')
     xlim([0 101])
     ylim([1e-10 1])
     %grid on
     ylabel('|H_y| (A/m)')
     %legend(['f=1 kHz Analytic'],['f=100 kHz Analytic'],['f=10 MHz Analytic'], ['f=1kHz FEM'], ['f=100 kHz FEM'], ['f=10 MHz FEM'])
     legend(['f=1 kHz Numerical integration'],['f=100 kHz Numerical integration'],['f=10 MHz Numerical integration'], ['f=1kHz FEM computation'], ['f=100 kHz FEM computation'], ['f=10 MHz FEM computation'])


     figure(12)
    %Hz
    semilogy(Rz(i,stepMAT),Hz1kMAT_y2(i,stepMAT),'sb')
    hold on
    semilogy(Rz(i,stepMAT),Hz100kMAT_y2(i,stepMAT),'dk')
    semilogy(Rz(i,stepMAT),Hz10MMAT_y2(i,stepMAT),'or')
    x_comsol=comsol_sol1kz(:,1);
    y_comsol=comsol_sol1kz(:,2);
    z_comsol=comsol_sol1kz(:,3);
    R_comsol=sqrt(x_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+y_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(z_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol)-he).^2);
    semilogy(R_comsol(rangeRcomsol),abs(Hz_comsol((i-1)*rx_rho_comsol+first:i*rx_rho_comsol,1)),'-b','Linewidth',1.5)
    semilogy(R_comsol(rangeRcomsol),abs(Hz_comsol((i-1)*rx_rho_comsol+first:i*rx_rho_comsol,2)),':k','Linewidth',1.5)
    semilogy(R_comsol(rangeRcomsol),abs(Hz_comsol((i-1)*rx_rho_comsol+first:i*rx_rho_comsol,3)),'--r','Linewidth',1.5)

     set(gca,'FontSize',12,'Fontname','arial')
     xlabel('R (m)')
     xlim([0 101])
     ylim([1e-10 1])
     %grid on
     ylabel('|H_z| (A/m)')
     %legend(['f=1 kHz Analytic'],['f=100 kHz Analytic'],['f=10 MHz Analytic'], ['f=1kHz FEM'], ['f=100 kHz FEM'], ['f=10 MHz FEM'])
     legend(['f=1 kHz Numerical integration'],['f=100 kHz Numerical integration'],['f=10 MHz Numerical integration'], ['f=1kHz FEM computation'], ['f=100 kHz FEM computation'], ['f=10 MHz FEM computation'])



 i=2; %for he=-7
    R_comsol=sqrt(x_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+y_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+ (z_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol)-he).^2);
    for id=1:length(px)
        indice=find(R_comsol==R(i,id));
        Hx_1k_com(id)=abs(Hx_comsol((i-1)*rx_rho_comsol+indice,1));
        Hx_100k_com(id)=abs(Hx_comsol((i-1)*rx_rho_comsol+indice,2));
        Hx_10M_com(id)=abs(Hx_comsol((i-1)*rx_rho_comsol+indice,3));
        Hy_1k_com(id)=abs(Hy_comsol((i-1)*rx_rho_comsol+indice,1));
        Hy_100k_com(id)=abs(Hy_comsol((i-1)*rx_rho_comsol+indice,2));
        Hy_10M_com(id)=abs(Hy_comsol((i-1)*rx_rho_comsol+indice,3));
        Hz_1k_com(id)=abs(Hz_comsol((i-1)*rx_rho_comsol+indice,1));
        Hz_100k_com(id)=abs(Hz_comsol((i-1)*rx_rho_comsol+indice,2));
        Hz_10M_com(id)=abs(Hz_comsol((i-1)*rx_rho_comsol+indice,3));
    end
    Hx_1k_mat=Hx1kMAT_xy(i,:);
    Hx_100k_mat=Hx100kMAT_xy(i,:);
    Hx_10M_mat=Hx10MMAT_xy(i,:);
    Hy_1k_mat=Hy1kMAT_y(i,:);
    Hy_100k_mat=Hy100kMAT_y(i,:);
    Hy_10M_mat=Hy10MMAT_y(i,:);
    Hz_1k_mat=Hz1kMAT_y2(i,:);
    Hz_100k_mat=Hz100kMAT_y2(i,:);
    Hz_10M_mat=Hz10MMAT_y2(i,:);
    
    
    
    R(i,:)
    relative_error_1k_Hx=(Hx_1k_mat-Hx_1k_com)./Hx_1k_mat
    relative_error_100k_Hx=(Hx_100k_mat-Hx_100k_com)./Hx_100k_mat
    relative_error_10M_Hx=(Hx_10M_mat-Hx_10M_com)./Hx_10M_mat
    relative_error_1k_Hy=(Hy_1k_mat-Hy_1k_com)./Hy_1k_mat
    relative_error_100k_Hy=(Hy_100k_mat-Hy_100k_com)./Hy_100k_mat
    relative_error_10M_Hy=(Hy_10M_mat-Hy_10M_com)./Hy_10M_mat
    relative_error_1k_Hz=(Hz_1k_mat-Hz_1k_com)./Hz_1k_mat
    relative_error_100k_Hz=(Hz_100k_mat-Hz_100k_com)./Hz_100k_mat
    relative_error_10M_Hz=(Hz_10M_mat-Hz_10M_com)./Hz_10M_mat
    
    
    range=[2:length(rho)];
    NRMSE_1k_Hx=1-(norm(Hx_1k_mat(range)-Hx_1k_com(range))/norm(Hx_1k_mat(range)-mean(Hx_1k_mat(range))));
    NRMSE_100k_Hx=1-(norm(Hx_100k_mat(range)-Hx_100k_com(range))/norm(Hx_100k_mat(range)-mean(Hx_100k_mat(range))));
    NRMSE_10M_Hx=1-(norm(Hx_10M_mat(range)-Hx_10M_com(range))/norm(Hx_10M_mat(range)-mean(Hx_10M_mat(range))));
    NRMSE_1k_Hy=1-(norm(Hy_1k_mat(range)-Hy_1k_com(range))/norm(Hz_1k_mat(range)-mean(Hz_1k_mat(range))));
    NRMSE_100k_Hy=1-(norm(Hy_100k_mat(range)-Hy_100k_com(range))/norm(Hz_100k_mat(range)-mean(Hz_100k_mat(range))));
    NRMSE_10M_Hy=1-(norm(Hy_10M_mat(range)-Hy_10M_com(range))/norm(Hz_10M_mat(range)-mean(Hz_10M_mat(range))));
    NRMSE_1k_Hz=1-(norm(Hz_1k_mat(range)-Hz_1k_com(range))/norm(Hz_1k_mat(range)-mean(Hz_1k_mat(range))));
    NRMSE_100k_Hz=1-(norm(Hz_100k_mat(range)-Hz_100k_com(range))/norm(Hz_100k_mat(range)-mean(Hz_100k_mat(range))));
    NRMSE_10M_Hz=1-(norm(Hz_10M_mat(range)-Hz_10M_com(range))/norm(Hz_10M_mat(range)-mean(Hz_10M_mat(range))));
    
    
    figure(13)
    semilogy(rho,Hx_1k_mat,'-b+')
    hold on
    semilogy(rho,Hx_1k_com,'-r*')
    semilogy(rho,Hx_100k_mat,'-b+')
    semilogy(rho,Hx_100k_com,'-r*')
    semilogy(rho,Hx_10M_mat,'-b+')
    semilogy(rho,Hx_10M_com,'-r*')
    xlim([0 101])
    ylim([1e-10 1])
    xlabel ('Hx')
  
    figure(14)
    semilogy(rho,Hy_1k_mat,'-b+')
    hold on
    semilogy(rho,Hy_1k_com,'-r*')
    semilogy(rho,Hy_100k_mat,'-b+')
    semilogy(rho,Hy_100k_com,'-r*')
    semilogy(rho,Hy_10M_mat,'-b+')
    semilogy(rho,Hy_10M_com,'-r*')
    xlim([0 101])
    ylim([1e-10 1])
    xlabel ('Hy')
  
    figure(15)
    semilogy(rho,Hz_1k_mat,'-b+')
    hold on
    semilogy(rho,Hz_1k_com,'-r*')
    semilogy(rho,Hz_100k_mat,'-b+')
    semilogy(rho,Hz_100k_com,'-r*')
    semilogy(rho,Hz_10M_mat,'-b+')
    semilogy(rho,Hz_10M_com,'-r*')
    xlim([0 101])
    ylim([1e-10 1])
    xlabel ('Hz')
  
     
    figure(16)
    plot(R(i,:),100*abs(relative_error_1k_Hx),'-b+')
    hold on
    semilogy(R(i,:),100*abs(relative_error_100k_Hx),'-k+')
    semilogy(R(i,:),100*abs(relative_error_10M_Hx),'-r+')
    legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])
    ylabel('Relative error (%) Hx')

    figure(17)
    plot(R(i,:),100*abs(relative_error_1k_Hy),'-b+')
    hold on
    semilogy(R(i,:),100*abs(relative_error_100k_Hy),'-k+')
    semilogy(R(i,:),100*abs(relative_error_10M_Hy),'-r+')
    legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])
    ylabel('Relative error (%) Hy')

    figure(18)
    plot(R(i,:),100*abs(relative_error_1k_Hz),'-b+')
    hold on
    semilogy(R(i,:),100*abs(relative_error_100k_Hz),'-k+')
    semilogy(R(i,:),100*abs(relative_error_10M_Hz),'-r+') 
    legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])
    ylabel('Relative error (%) Hz')
    
    fprintf('Mean relative error 1kHz Hx: %d\n',mean(relative_error_1k_Hx(2:length(rho))));
    fprintf('Mean relative error 100kHz Hx: %d\n',mean(relative_error_100k_Hx(2:length(rho))));
    fprintf('Mean relative error 10MHz Hx: %d\n',mean(relative_error_10M_Hx(2:length(rho))));
    fprintf('Mean relative error 1kHz Hy: %d\n',mean(relative_error_1k_Hy(2:length(rho))));
    fprintf('Mean relative error 100kHz Hy: %d\n',mean(relative_error_100k_Hy(2:length(rho))));
    fprintf('Mean relative error 10MHz Hy: %d\n',mean(relative_error_10M_Hy(2:length(rho))));
    fprintf('Mean relative error 1kHz Hz: %d\n',mean(relative_error_1k_Hz(2:length(rho))));
    fprintf('Mean relative error 100kHz Hz: %d\n',mean(relative_error_100k_Hz(2:length(rho))));
    fprintf('Mean relative error 10MHz Hz: %d\n',mean(relative_error_10M_Hz(2:length(rho))));
    
    
    fprintf('Mean relative absolute error 1kHz Hx: %d\n',mean(abs(relative_error_1k_Hx(2:length(rho)))));
    fprintf('Mean relative absolute error 100kHz Hx: %d\n',mean(abs(relative_error_100k_Hx(2:length(rho)))));
    fprintf('Mean relative absolute error 10MHz Hx: %d\n',mean(abs(relative_error_10M_Hx(2:length(rho)))));
    fprintf('Mean relative absolute error 1kHz Hy: %d\n',mean(abs(relative_error_1k_Hy(2:length(rho)))));
    fprintf('Mean relative absolute error 100kHz Hy: %d\n',mean(abs(relative_error_100k_Hy(2:length(rho)))));
    fprintf('Mean relative absolute error 10MHz Hy: %d\n',mean(abs(relative_error_10M_Hy(2:length(rho)))));
    fprintf('Mean relative absolute error 1kHz Hz: %d\n',mean(abs(relative_error_1k_Hz(2:length(rho)))));
    fprintf('Mean relative absolute error 100kHz Hz: %d\n',mean(abs(relative_error_100k_Hz(2:length(rho)))));
    fprintf('Mean relative absolute error 10MHz Hz: %d\n',mean(abs(relative_error_10M_Hz(2:length(rho)))));
    
    fprintf('Mean relative absolute error 1kHz Hx: %f ppm\n',1e7*mean(abs(relative_error_1k_Hx(2:length(rho)))));
    fprintf('Mean relative absolute error 100kHz Hx: %f ppm\n',1e7*mean(abs(relative_error_100k_Hx(2:length(rho)))));
    fprintf('Mean relative absolute error 10MHz Hx: %f  ppm\n',1e7*mean(abs(relative_error_10M_Hx(2:length(rho)))));
    fprintf('Mean relative absolute error 1kHz Hy: %f ppm\n',1e7*mean(abs(relative_error_1k_Hy(2:length(rho)))));
    fprintf('Mean relative absolute error 100kHz Hy: %f ppm\n',1e7*mean(abs(relative_error_100k_Hy(2:length(rho)))));
    fprintf('Mean relative absolute error 10MHz Hy: %f  ppm\n',1e7*mean(abs(relative_error_10M_Hy(2:length(rho)))));
    fprintf('Mean relative absolute error 1kHz Hz: %f  ppm\n',1e7*mean(abs(relative_error_1k_Hz(2:length(rho)))));
    fprintf('Mean relative absolute error 100kHz Hz: %f ppm\n',1e7*mean(abs(relative_error_100k_Hz(2:length(rho)))));
    fprintf('Mean relative absolute error 10MHz Hz: %f ppm\n',1e7*mean(abs(relative_error_10M_Hz(2:length(rho)))));
    
    fprintf('NRMSE 1kHz Hx: %f \n',NRMSE_1k_Hx)
    fprintf('NRMSE 100kHz Hx: %f \n',NRMSE_100k_Hx)
    fprintf('NRMSE 10MHz Hx: %f \n',NRMSE_10M_Hx)
    fprintf('NRMSE 1kHz Hy: %f \n',NRMSE_1k_Hy)
    fprintf('NRMSE 100kHz Hy: %f \n',NRMSE_100k_Hy)
    fprintf('NRMSE 10MHz Hy: %f \n',NRMSE_10M_Hy)
    fprintf('NRMSE 1kHz Hz: %f \n',NRMSE_1k_Hz)
    fprintf('NRMSE 100kHz Hz: %f \n',NRMSE_100k_Hz)
    fprintf('NRMSE 10MHz Hz: %f \n',NRMSE_10M_Hz)
