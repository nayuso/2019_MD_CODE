%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VALIDATION OF THE ANALYTICAL EXPRESSIONS DUE TO A VMD
% IN A THREE-LAYERED REGION BY MEANS OF FEM SIMULATIONS
%
% Covered range: From near field to far field
%
%
% validation_VMDlayer.m
%
% Author: Dr. Natalia Ayuso
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



close all
clear all
clc

global mu0 e0 c

e0=8.85e-12;
mu0=4*pi*1e-7;
c=1/sqrt(e0*mu0);

% Emitter
m0=1;
he=-7;

f=[1e3 100e3 1e7]
w = 2*pi*f;
a=1e-2;

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
rho=[0:1:100];
py=rho;
px=zeros(1,length(py));
rx_hr=length(hr_vector);
rx_rho=length(rho);

fprintf('Emitter placed in region 2 \n');
load( '../../MATLAB_FEM/layered_soil2_MATLAB_fine');
 
fprintf('Field computation time: %3.0f min %2.0f s\n',minutes,segs);

for i=1:rx_hr
    R(i,:)=sqrt(rho.^2+(hr_vector(i)-he)^2);
    Hr_1k_MAT(i,:)= H_MAT1k(1,((i-1)*rx_rho)+1:rx_rho*i);
    Hz_1k_MAT(i,:)= H_MAT1k(3,((i-1)*rx_rho)+1:rx_rho*i);
    Hr_100k_MAT(i,:)= H_MAT100k(1,((i-1)*rx_rho)+1:rx_rho*i);
    Hz_100k_MAT(i,:)= H_MAT100k(3,((i-1)*rx_rho)+1:rx_rho*i);
    Hr_10M_MAT(i,:)= H_MAT10M(1,((i-1)*rx_rho)+1:rx_rho*i);
    Hz_10M_MAT(i,:)= H_MAT10M(3,((i-1)*rx_rho)+1:rx_rho*i);
end

Hr1kMAT=abs(Hr_1k_MAT);
Hr100kMAT=abs(Hr_100k_MAT);
Hr10MMAT=abs(Hr_10M_MAT);

Hz1kMAT=abs(Hz_1k_MAT);
Hz100kMAT=abs(Hz_100k_MAT);
Hz10MMAT=abs(Hz_10M_MAT);
     
style1k=['-+b -.+b--+b']; %A line style for any receiver heigth and frequency
style100k=['-*g -.*g--*g']; %A line style for any receiver heigth and frequency
style10M=['-or -.or--or']; %A line style for any receiver heigth and frequency

fprintf('Solid line for hr=1, dotted for hr=-2 and dashed for hr=-7\n')

figure(1)
for i=1:rx_hr
    semilogy(R(i,:),Hr1kMAT(i,:),style1k([(i-1)*4+1:i*4]),'LineWidth',2)
    hold on
    semilogy(R(i,:),Hr100kMAT(i,:),style100k([(i-1)*4+1:i*4]),'LineWidth',2)
    semilogy(R(i,:),Hr10MMAT(i,:),style10M([(i-1)*4+1:i*4]),'LineWidth',2)
end
set(gca,'FontSize',12,'Fontname','arial')
xlabel('R')
grid on
ylabel('|H_r| (A/m)')
title('MATLAB SOLUTION')
legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])

figure(2)
for i=1:rx_hr
    semilogy(R(i,:),Hz1kMAT(i,:),style1k([(i-1)*4+1:i*4]),'LineWidth',2)
    hold on
    semilogy(R(i,:),Hz100kMAT(i,:),style100k([(i-1)*4+1:i*4]),'LineWidth',2)
    semilogy(R(i,:),Hz10MMAT(i,:),style10M([(i-1)*4+1:i*4]),'LineWidth',2)
end
set(gca,'FontSize',12,'Fontname','arial')
xlabel('R')
grid on
ylabel('|H_z|')
title('MATLAB SOLUTION')
legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])

%IMPORT COMSOL DATA. Don't forguet to write down the file extention
fprintf('Emitter placed in region 2 \n');
comsol_sol_1k=load( '../../MATLAB_FEM/VMD_Tx_soil2_COMSOL_1k_mesh4_tol4.txt');
comsol_sol_100k=load( '../../MATLAB_FEM/VMD_Tx_soil2_COMSOL_mesh_8elementspoints.txt');
comsol_sol_10M=load( '../../MATLAB_FEM/VMD_Tx_soil2_COMSOL_10M_mesh13_tol3.txt');
    

r_comsol=comsol_sol(:,1);
z_comsol=comsol_sol(:,2);
    
Hr_comsol=[comsol_sol(:,3) comsol_sol(:,5) comsol_sol(:,7)];
Hz_comsol=[comsol_sol(:,4) comsol_sol(:,6) comsol_sol(:,8)];
    %Hr_comsol=[comsol_sol(:,3) comsol_sol_high(:,3) comsol_sol_high(:,5)];
    %Hz_comsol=[comsol_sol(:,4) comsol_sol_high(:,4) comsol_sol_high(:,6)];
    %Hr_comsol=[comsol_sol_1k(:,3) comsol_sol_100k(:,3) comsol_sol_10M(:,3)];
    %Hz_comsol=[comsol_sol_1k(:,4) comsol_sol_100k(:,4) comsol_sol_10M(:,4)];
       
    py_comsol=[0:100];
    px_comsol=zeros(1,length(py_comsol));
    rho_comsol=sqrt(px_comsol.^2+py_comsol.^2);
    rx_rho_comsol=length(rho_comsol);
    hr_vector_comsol=[1 -2 -7];
    %hr_vector_comsol=[1];
    rx_hr_comsol=length(hr_vector_comsol);

    style_comsol=['- -.--'];

    figure(3)
    for i=1:rx_hr_comsol
        R_comsol=sqrt(r_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(hr_vector_comsol(i)-he)^2);
        semilogy(R_comsol,abs(Hr_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,:)),style_comsol([(i-1)*2+1:i*2]))
        hold on
    end
    set(gca,'FontSize',12,'Fontname','arial')
    xlabel('R')
    grid on
    ylabel('|H_r|')
    title('COMSOL SOLUTION')
    legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])

    figure(4)
    for i=1:rx_hr_comsol
        R_comsol=sqrt(r_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(hr_vector_comsol(i)-he)^2);
        semilogy(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,:)),style_comsol([(i-1)*2+1:i*2]))
        hold on
    end
    set(gca,'FontSize',12,'Fontname','arial')
    xlabel('R')
    grid on
    ylabel('|H_Z|')
    title('COMSOL SOLUTION')
    legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])

    figure(5)
    for i=1:rx_hr
        loglog(R(i,:),Hr1kMAT(i,:),'+b')
        hold on
        loglog(R(i,:),Hr100kMAT(i,:),'sg')
        loglog(R(i,:),Hr10MMAT(i,:),'or')
    end
    for i=1:rx_hr_comsol
        R_comsol=sqrt(r_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(hr_vector_comsol(i)-he)^2);
        loglog(R_comsol,abs(Hr_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,:)),style_comsol([(i-1)*2+1:i*2]))
    end
    set(gca,'FontSize',12,'Fontname','arial')
    xlabel('R')
    xlim([0 100])
    grid on
    ylabel('|H_r| (A/m)')
    title('MATLAB COMSOL SOLUTION')
    legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])

        figure(6)
    for i=1:rx_hr
        loglog(R(i,:),Hz1kMAT(i,:),'+b')
        hold on
        loglog(R(i,:),Hz100kMAT(i,:),'sg')
        loglog(R(i,:),Hz10MMAT(i,:),'or')
    end
    for i=1:rx_hr_comsol
        R_comsol=sqrt(r_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(hr_vector_comsol(i)-he)^2);
        loglog(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,:)),style_comsol([(i-1)*2+1:i*2]))
    end
    set(gca,'FontSize',12,'Fontname','arial')
    xlabel('R')
    xlim([0 100])
    grid on
    ylabel('|H_Z| (A/m)')
    title('MATLAB COMSOL SOLUTION')
    legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])
else
    return
end

order1k=[4 7 10];
order100k=[4 7 10];
order10M=[4 7 10];
figure(7)
    for i=1:rx_hr
        loglog(R(i,:),Hr1kMAT(i,:),'vb')
        hold on
        loglog(R(i,:),Hr100kMAT(i,:),'dk')
        loglog(R(i,:),Hr10MMAT(i,:),'or')
        R_comsol=sqrt(r_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(hr_vector_comsol(i)-he)^2);
        loglog(R_comsol,abs(Hr_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,1)),'-b','Linewidth',1.5)
         text(R(i,order1k(i)),Hr1kMAT(i,order1k(i)),['z =  ',num2str(hr_vector(i)),'\rightarrow'],...
 'VerticalAlignment','middle',...
 'HorizontalAlignment','right',...
  'FontSize',14)
%'Backgroundcolor','white',...
  %'EdgeColor','black',...
        loglog(R_comsol,abs(Hr_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,2)),':k','Linewidth',1.5)
        text(R(i,order100k(i)),Hr100kMAT(i,order100k(i)),['z =  ',num2str(hr_vector(i)),'\rightarrow'],...
 'VerticalAlignment','middle',...
 'HorizontalAlignment','right',...
 'FontSize',14)
        loglog(R_comsol,abs(Hr_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,3)),'--r','Linewidth',1.5)      
        text(R(i,order10M(i)),Hr10MMAT(i,order10M(i)),['z =  ',num2str(hr_vector(i)),'\rightarrow'],...
 'VerticalAlignment','middle',...
 'HorizontalAlignment','right',...
 'FontSize',14)
    end
    
    set(gca,'FontSize',12,'Fontname','arial')
    xlabel('R (m)')
    xlim([1 100])
    ylim([1e-10 1])
    %grid on
    ylabel('|H_r| (A/m)')
    %legend(['f=1 kHz MATLAB'],['f=100 kHz MATLAB'],['f=10 MHz MATLAB'], ['f=1kHz COMSOL'], ['f=100 kHz COMSOL'], ['f=10 MHz COMSOL'])
    legend(['f=1 kHz Analytic'],['f=100 kHz Analytic'],['f=10 MHz Analytic'], ['f=1kHz FEM'], ['f=100 kHz FEM'], ['f=10 MHz FEM'])
    
    order1k=[4 7 10];
order100k=[4 7 10];
order10M=[4 7 10];
figure(8)
    for i=1:rx_hr
      loglog(R(i,:),Hz1kMAT(i,:),'vb')
        hold on
        loglog(R(i,:),Hz100kMAT(i,:),'dk')
        loglog(R(i,:),Hz10MMAT(i,:),'or')
        R_comsol=sqrt(r_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(hr_vector_comsol(i)-he)^2);
        loglog(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,1)),'-b','Linewidth',1.5)
         text(R(i,order1k(i)),Hz1kMAT(i,order1k(i)),['z =  ',num2str(hr_vector(i)),'\rightarrow'],...
 'VerticalAlignment','middle',...
 'HorizontalAlignment','right',...
  'FontSize',14)
%'Backgroundcolor','white',...
  %'EdgeColor','black',...
        loglog(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,2)),':k','Linewidth',1.5)
        text(R(i,order100k(i)),Hz100kMAT(i,order100k(i)),['z =  ',num2str(hr_vector(i)),'\rightarrow'],...
 'VerticalAlignment','middle',...
 'HorizontalAlignment','right',...
 'FontSize',14)
        loglog(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,3)),'--r','Linewidth',1.5)      
        text(R(i,order10M(i)),Hz10MMAT(i,order10M(i)),['z =  ',num2str(hr_vector(i)),'\rightarrow'],...
 'VerticalAlignment','middle',...
 'HorizontalAlignment','right',...
 'FontSize',14)
    end
    
    set(gca,'FontSize',12,'Fontname','arial')
    xlabel('R (m)')
    %xlim([1 100])
    xlim([1 101])
    ylim([1e-10 1])
    %grid on
    ylabel('|H_z| (A/m)')
    %legend(['f=1 kHz MATLAB'],['f=100 kHz MATLAB'],['f=10 MHz MATLAB'], ['f=1kHz COMSOL'], ['f=100 kHz COMSOL'], ['f=10 MHz COMSOL'])
    legend(['f=1 kHz Numerical integration'],['f=100 kHz Numerical integration'],['f=10 MHz Numerical integration'], ['f=1kHz FEM computation'], ['f=100 kHz FEM computation'], ['f=10 MHz FEM computation'])
    
    %Tx=1 Rx=-2
    figure(9)
    i=2
    semilogy(R(i,:),Hr1kMAT(i,:),'vb')
    hold on
    semilogy(R(i,:),Hr100kMAT(i,:),'dk')
    semilogy(R(i,:),Hr10MMAT(i,:),'or')
    R_comsol=sqrt(r_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(hr_vector_comsol(i)-he)^2);
    semilogy(R_comsol(2:length(R_comsol)),abs(Hr_comsol((i-1)*rx_rho_comsol+2:i*rx_rho_comsol,1)),'-b','Linewidth',1.5)
    semilogy(R_comsol(2:length(R_comsol)),abs(Hr_comsol((i-1)*rx_rho_comsol+2:i*rx_rho_comsol,2)),':k','Linewidth',1.5)
    semilogy(R_comsol(2:length(R_comsol)),abs(Hr_comsol((i-1)*rx_rho_comsol+2:i*rx_rho_comsol,3)),'--r','Linewidth',1.5)
            
    set(gca,'FontSize',12,'Fontname','arial')
    xlabel('R (m)')
    %xlim([0 100])
    xlim([1 101])
    ylim([1e-10 1])
    %grid on
    ylabel('|H_r| (A/m)')
    %legend(['f=1 kHz MATLAB'],['f=100 kHz MATLAB'],['f=10 MHz MATLAB'], ['f=1kHz COMSOL'], ['f=100 kHz COMSOL'], ['f=10 MHz COMSOL'])
    %legend(['f=1 kHz Analytic'],['f=100 kHz Analytic'],['f=10 MHz Analytic'], ['f=1kHz FEM'], ['f=100 kHz FEM'], ['f=10 MHz FEM'])
    legend(['f=1 kHz Numerical integration'],['f=100 kHz Numerical integration'],['f=10 MHz Numerical integration'], ['f=1kHz FEM computation'], ['f=100 kHz FEM computation'], ['f=10 MHz FEM computation'])
    
    %Tx=1 Rx=-2
    figure(10)
    semilogy(R(i,:),Hz1kMAT(i,:),'vb')
    hold on
    semilogy(R(i,:),Hz100kMAT(i,:),'dk')
    semilogy(R(i,:),Hz10MMAT(i,:),'or')
    R_comsol=sqrt(r_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(hr_vector_comsol(i)-he)^2);
    semilogy(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,1)),'-b','Linewidth',1.5)
    semilogy(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,2)),':k','Linewidth',1.5)
    semilogy(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,3)),'--r','Linewidth',1.5)
            
    set(gca,'FontSize',12,'Fontname','arial')
    xlabel('R (m)')
    %xlim([0 100])
    xlim([1 101])
    ylim([1e-10 1])
    %grid on
    ylabel('|H_z| (A/m)')
    %legend(['f=1 kHz MATLAB'],['f=100 kHz MATLAB'],['f=10 MHz MATLAB'], ['f=1kHz COMSOL'], ['f=100 kHz COMSOL'], ['f=10 MHz COMSOL'])
    %legend(['f=1 kHz Analytic'],['f=100 kHz Analytic'],['f=10 MHz Analytic'], ['f=1kHz FEM'], ['f=100 kHz FEM'], ['f=10 MHz FEM'])
    legend(['f=1 kHz Numerical integration'],['f=100 kHz Numerical integration'],['f=10 MHz Numerical integration'], ['f=1kHz FEM computation'], ['f=100 kHz FEM computation'], ['f=10 MHz FEM computation'])
    
    stepMAT=[1:2:101];
    %Rx=1(i=1)
    figure(11)
    i=1
    semilogy(R(i,stepMAT),Hr1kMAT(i,stepMAT),'vb')
    hold on
    semilogy(R(i,stepMAT),Hr100kMAT(i,stepMAT),'dk')
    semilogy(R(i,stepMAT),Hr10MMAT(i,stepMAT),'or')
    R_comsol=sqrt(r_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(hr_vector_comsol(i)-he)^2);
    semilogy(R_comsol(2:length(R_comsol)),abs(Hr_comsol((i-1)*rx_rho_comsol+2:i*rx_rho_comsol,1)),'-b','Linewidth',1.5)
    semilogy(R_comsol(2:length(R_comsol)),abs(Hr_comsol((i-1)*rx_rho_comsol+2:i*rx_rho_comsol,2)),':k','Linewidth',1.5)
    semilogy(R_comsol(2:length(R_comsol)),abs(Hr_comsol((i-1)*rx_rho_comsol+2:i*rx_rho_comsol,3)),'--r','Linewidth',1.5)
            
    set(gca,'FontSize',12,'Fontname','arial')
    xlabel('R (m)')
    %xlim([0 100])
    xlim([0 101])
    ylim([1e-10 1])
    %grid on
    ylabel('|H_r| (A/m)')
    %legend(['f=1 kHz MATLAB'],['f=100 kHz MATLAB'],['f=10 MHz MATLAB'], ['f=1kHz COMSOL'], ['f=100 kHz COMSOL'], ['f=10 MHz COMSOL'])
    %legend(['f=1 kHz Analytic'],['f=100 kHz Analytic'],['f=10 MHz Analytic'], ['f=1kHz FEM'], ['f=100 kHz FEM'], ['f=10 MHz FEM'])
    legend(['f=1 kHz Numerical integration'],['f=100 kHz Numerical integration'],['f=10 MHz Numerical integration'], ['f=1kHz FEM computation'], ['f=100 kHz FEM computation'], ['f=10 MHz FEM computation'])
    
    %Rx=-2 (i=1)
    figure(12)
    semilogy(R(i,stepMAT),Hz1kMAT(i,stepMAT),'vb')
    hold on
    semilogy(R(i,stepMAT),Hz100kMAT(i,stepMAT),'dk')
    semilogy(R(i,stepMAT),Hz10MMAT(i,stepMAT),'or')
    R_comsol=sqrt(r_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(hr_vector_comsol(i)-he)^2);
    semilogy(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,1)),'-b','Linewidth',1.5)
    semilogy(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,2)),':k','Linewidth',1.5)
    semilogy(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,3)),'--r','Linewidth',1.5)
            
    set(gca,'FontSize',12,'Fontname','arial')
    xlabel('R (m)')
    %xlim([0 100])
    xlim([0 101])
    ylim([1e-10 1])
    %grid on
    ylabel('|H_z| (A/m)')
    %legend(['f=1 kHz MATLAB'],['f=100 kHz MATLAB'],['f=10 MHz MATLAB'], ['f=1kHz COMSOL'], ['f=100 kHz COMSOL'], ['f=10 MHz COMSOL'])
    %legend(['f=1 kHz Analytic'],['f=100 kHz Analytic'],['f=10 MHz Analytic'], ['f=1kHz FEM'], ['f=100 kHz FEM'], ['f=10 MHz FEM'])
    legend(['f=1 kHz Numerical integration'],['f=100 kHz Numerical integration'],['f=10 MHz Numerical integration'], ['f=1kHz FEM computation'], ['f=100 kHz FEM computation'], ['f=10 MHz FEM computation'])
    
    %Rx=-2 (i=2)
    figure(13)
    i=2
    semilogy(R(i,stepMAT),Hr1kMAT(i,stepMAT),'vb')
    hold on
    semilogy(R(i,stepMAT),Hr100kMAT(i,stepMAT),'dk')
    semilogy(R(i,stepMAT),Hr10MMAT(i,stepMAT),'or')
    R_comsol=sqrt(r_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(hr_vector_comsol(i)-he)^2);
    semilogy(R_comsol(2:length(R_comsol)),abs(Hr_comsol((i-1)*rx_rho_comsol+2:i*rx_rho_comsol,1)),'-b','Linewidth',1.5)
    semilogy(R_comsol(2:length(R_comsol)),abs(Hr_comsol((i-1)*rx_rho_comsol+2:i*rx_rho_comsol,2)),':k','Linewidth',1.5)
    semilogy(R_comsol(2:length(R_comsol)),abs(Hr_comsol((i-1)*rx_rho_comsol+2:i*rx_rho_comsol,3)),'--r','Linewidth',1.5)
            
    set(gca,'FontSize',12,'Fontname','arial')
    xlabel('R (m)')
    %xlim([0 100])
    xlim([0 101])
    ylim([1e-10 1])
    %grid on
    ylabel('|H_r| (A/m)')
    %legend(['f=1 kHz MATLAB'],['f=100 kHz MATLAB'],['f=10 MHz MATLAB'], ['f=1kHz COMSOL'], ['f=100 kHz COMSOL'], ['f=10 MHz COMSOL'])
    %legend(['f=1 kHz Analytic'],['f=100 kHz Analytic'],['f=10 MHz Analytic'], ['f=1kHz FEM'], ['f=100 kHz FEM'], ['f=10 MHz FEM'])
    legend(['f=1 kHz Numerical integration'],['f=100 kHz Numerical integration'],['f=10 MHz Numerical integration'], ['f=1kHz FEM computation'], ['f=100 kHz FEM computation'], ['f=10 MHz FEM computation'])
    
    %Rx=-2
    figure(14)
    semilogy(R(i,stepMAT),Hz1kMAT(i,stepMAT),'vb')
    hold on
    semilogy(R(i,stepMAT),Hz100kMAT(i,stepMAT),'dk')
    semilogy(R(i,stepMAT),Hz10MMAT(i,stepMAT),'or')
    R_comsol=sqrt(r_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(hr_vector_comsol(i)-he)^2);
    semilogy(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,1)),'-b','Linewidth',1.5)
    semilogy(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,2)),':k','Linewidth',1.5)
    semilogy(R_comsol,abs(Hz_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol,3)),'--r','Linewidth',1.5)
            
    set(gca,'FontSize',12,'Fontname','arial')
    xlabel('R (m)')
    %xlim([0 100])
    xlim([0 101])
    ylim([1e-10 1])
    %grid on
    ylabel('|H_z| (A/m)')
    %legend(['f=1 kHz MATLAB'],['f=100 kHz MATLAB'],['f=10 MHz MATLAB'], ['f=1kHz COMSOL'], ['f=100 kHz COMSOL'], ['f=10 MHz COMSOL'])
    %legend(['f=1 kHz Analytic'],['f=100 kHz Analytic'],['f=10 MHz Analytic'], ['f=1kHz FEM'], ['f=100 kHz FEM'], ['f=10 MHz FEM'])
    legend(['f=1 kHz Numerical integration'],['f=100 kHz Numerical integration'],['f=10 MHz Numerical integration'], ['f=1kHz FEM computation'], ['f=100 kHz FEM computation'], ['f=10 MHz FEM computation'])
    
   switch he
       case -2
        i=1;
       case 1
            i=2;
       case -7
           i=2;
  end
    R_comsol=sqrt(r_comsol((i-1)*rx_rho_comsol+1:i*rx_rho_comsol).^2+(hr_vector_comsol(i)-he)^2);
    for id=1:length(y)
        indice=find(R_comsol==R(i,id));
        Hr_1k_com(id)=abs(Hr_comsol((i-1)*rx_rho_comsol+indice,1));
        Hr_100k_com(id)=abs(Hr_comsol((i-1)*rx_rho_comsol+indice,2));
        Hr_10M_com(id)=abs(Hr_comsol((i-1)*rx_rho_comsol+indice,3));
        Hz_1k_com(id)=abs(Hz_comsol((i-1)*rx_rho_comsol+indice,1));
        Hz_100k_com(id)=abs(Hz_comsol((i-1)*rx_rho_comsol+indice,2));
        Hz_10M_com(id)=abs(Hz_comsol((i-1)*rx_rho_comsol+indice,3));
    end
    Hr_1k_mat=Hr1kMAT(i,:);
    Hr_100k_mat=Hr100kMAT(i,:);
    Hr_10M_mat=Hr10MMAT(i,:);
    Hz_1k_mat=Hz1kMAT(i,:);
    Hz_100k_mat=Hz100kMAT(i,:);
    Hz_10M_mat=Hz10MMAT(i,:);
    
    R(i,:)
    relative_error_1k_Hr=(Hr_1k_mat-Hr_1k_com)./Hr_1k_com
    relative_error_100k_Hr=(Hr_100k_mat-Hr_100k_com)./Hr_100k_com
    relative_error_10M_Hr=(Hr_10M_mat-Hr_10M_com)./Hr_10M_com
    relative_error_1k_Hz=(Hz_1k_mat-Hz_1k_com)./Hz_1k_com
    relative_error_100k_Hz=(Hz_100k_mat-Hz_100k_com)./Hz_100k_com
    relative_error_10M_Hz=(Hz_10M_mat-Hz_10M_com)./Hz_10M_com
    
    rmse_1k_Hr=(1/(length(rho)-1))*sum(Hr_1k_mat(2:length(rho))-Hr_1k_com(2:length(rho))).^2
    rmse_100k_Hr=(1/(length(rho)-1))*sum(Hr_100k_mat(2:length(rho))-Hr_100k_com(2:length(rho))).^2
    rmse_10M_Hr=(1/(length(rho)-1))*sum(Hr_10M_mat(2:length(rho))-Hr_10M_com(2:length(rho))).^2
    rmse_1k_Hz=1/length(rho)*sum((Hz_1k_mat-Hz_1k_com).^2)
    rmse_100k_Hz=1/length(rho)*sum((Hz_100k_mat-Hz_100k_com).^2)
    rmse_10M_Hz=1/length(rho)*sum((Hz_10M_mat-Hz_10M_com).^2)
    
    NRMSE_1k_Hr=1-(norm(Hr_1k_mat-Hr_1k_com)/norm(Hr_1k_mat-mean(Hr_1k_mat)));
    NRMSE_100k_Hr=1-(norm(Hr_100k_mat-Hr_100k_com)/norm(Hr_100k_mat-mean(Hr_100k_mat)));
    NRMSE_10M_Hr=1-(norm(Hr_10M_mat-Hr_10M_com)/norm(Hr_10M_mat-mean(Hr_10M_mat)));
    NRMSE_1k_Hz=1-(norm(Hz_1k_mat-Hz_1k_com)/norm(Hz_1k_mat-mean(Hz_1k_mat)));
    NRMSE_100k_Hz=1-(norm(Hz_100k_mat-Hz_100k_com)/norm(Hz_100k_mat-mean(Hz_100k_mat)));
    NRMSE_10M_Hz=1-(norm(Hz_10M_mat-Hz_10M_com)/norm(Hz_10M_mat-mean(Hz_10M_mat)));
          
    figure(15)
    semilogy(y,Hr_1k_mat,'-b+')
    hold on
    semilogy(y,Hr_1k_com,'-r*')
    semilogy(y,Hr_100k_mat,'-b+')
    semilogy(y,Hr_100k_com,'-r*')
    semilogy(y,Hr_10M_mat,'-b+')
    semilogy(y,Hr_10M_com,'-r*')
    xlim([0 101])
    ylim([1e-10 1])
    xlabel ('Hr')
  
    figure(16)
    semilogy(y,Hz_1k_mat,'-b+')
    hold on
    semilogy(y,Hz_1k_com,'-r*')
    semilogy(y,Hz_100k_mat,'-b+')
    semilogy(y,Hz_100k_com,'-r*')
    semilogy(y,Hz_10M_mat,'-b+')
    semilogy(y,Hz_10M_com,'-r*')
    xlim([0 101])
    ylim([1e-10 1])
    xlabel ('Hr')
     
    figure(17)
    plot(R(i,2:length(rho)),100*abs(relative_error_1k_Hr(2:length(rho))),'-b+')
    hold on
    semilogy(R(i,2:length(rho)),100*abs(relative_error_100k_Hr(2:length(rho))),'-k+')
    semilogy(R(i,2:length(rho)),100*abs(relative_error_10M_Hr(2:length(rho))),'-r+')
    legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])
    ylabel('Relative error (%) Hr')

    figure(18)
    plot(R(i,:),100*abs(relative_error_1k_Hz),'-b+')
    hold on
    semilogy(R(i,:),100*abs(relative_error_100k_Hz),'-k+')
    semilogy(R(i,:),100*abs(relative_error_10M_Hz),'-r+') 
    legend(['f=' num2str(f(1))],['f=' num2str(f(2))],['f=' num2str(f(3))])
    ylabel('Relative error (%) Hz')
    
    fprintf('Mean relative error 1kHz Hr: %d\n',mean(relative_error_1k_Hr(2:length(rho))));
    fprintf('Mean relative error 100kHz Hr: %d\n',mean(relative_error_100k_Hr(2:length(rho))));
    fprintf('Mean relative error 10MHz Hr: %d\n',mean(relative_error_10M_Hr(2:length(rho))));
    fprintf('Mean relative error 1kHz Hz: %d\n',mean(relative_error_1k_Hz(1:length(rho))));
    fprintf('Mean relative error 100kHz Hz: %d\n',mean(relative_error_100k_Hz(1:length(rho))));
    fprintf('Mean relative error 10MHz Hz: %d\n',mean(relative_error_10M_Hz(1:length(rho))));
    
    fprintf('Mean relative absolute error 1kHz Hr: %d\n',mean(abs(relative_error_1k_Hr(2:length(rho)))));
    fprintf('Mean relative absolute error 100kHz Hr: %d\n',mean(abs(relative_error_100k_Hr(2:length(rho)))));
    fprintf('Mean relative absolute error 10MHz Hr: %d\n',mean(abs(relative_error_10M_Hr(2:length(rho)))));
    fprintf('Mean relative absolute error 1kHz Hz: %d\n',mean(abs(relative_error_1k_Hz(1:length(rho)))));
    fprintf('Mean relative absolute error 100kHz Hz: %d\n',mean(abs(relative_error_100k_Hz(1:length(rho)))));
    fprintf('Mean relative absolute error 10MHz Hz: %d\n',mean(abs(relative_error_10M_Hz(1:length(rho)))));
    
    fprintf('Mean relative absolute error 1kHz Hr: %f ppm\n',1e7*mean(abs(relative_error_1k_Hr(2:length(rho)))));
    fprintf('Mean relative absolute error 100kHz Hr: %f ppm\n',1e7*mean(abs(relative_error_100k_Hr(2:length(rho)))));
    fprintf('Mean relative absolute error 10MHz Hr: %f  ppm\n',1e7*mean(abs(relative_error_10M_Hr(2:length(rho)))));
    fprintf('Mean relative absolute error 1kHz Hz: %f  ppm\n',1e7*mean(abs(relative_error_1k_Hz(1:length(rho)))));
    fprintf('Mean relative absolute error 100kHz Hz: %f ppm\n',1e7*mean(abs(relative_error_100k_Hz(1:length(rho)))));
    fprintf('Mean relative absolute error 10MHz Hz: %f ppm\n',1e7*mean(abs(relative_error_10M_Hz(1:length(rho)))));
    
    fprintf('RMSE 1kHz Hr: %d\n',rmse_1k_Hr);
    fprintf('RMSE 100kHz Hr: %d\n',rmse_100k_Hr);
    fprintf('RMSE 10MHz Hr: %d\n',rmse_10M_Hr);
    fprintf('RMSE 1kHz Hz: %d\n',rmse_1k_Hz);
    fprintf('RMSE 100kHz Hz: %d\n',rmse_100k_Hz);
    fprintf('RMSE 10MHz Hz: %d\n',rmse_10M_Hz);
    
    fprintf('RMSE 1kHz Hr: %f ppm\n',1e7*rmse_1k_Hr);
    fprintf('RMSE 100kHz Hr: %f ppm\n',1e7*rmse_100k_Hr);
    fprintf('RMSE 10MHz Hr: %f  ppm\n',1e7*rmse_10M_Hr);
    fprintf('RMSE 1kHz Hz: %f  ppm\n',1e7*rmse_1k_Hz);
    fprintf('RMSE 100kHz Hz: %f ppm\n',1e7*rmse_100k_Hz);
    fprintf('RMSE 10MHz Hz: %f ppm\n',1e7*rmse_10M_Hz);
    
    fprintf('NRMSE 1kHz Hr: %f \n',NRMSE_1k_Hr)
    fprintf('NRMSE 100kHz Hr: %f \n',NRMSE_100k_Hr)
    fprintf('NRMSE 10MHz Hr: %f \n',NRMSE_10M_Hr)
    fprintf('NRMSE 1kHz Hz: %f \n',NRMSE_1k_Hz)
    fprintf('NRMSE 100kHz Hz: %f \n',NRMSE_100k_Hz)
    fprintf('NRMSE 10MHz Hz: %f \n',NRMSE_10M_Hz)
