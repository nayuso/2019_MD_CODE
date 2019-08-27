%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Application: sensor nodes that wireless communicate under the ground
% The best operating frequency is still an open research aspect
% According to Silva and Moghaddam: an approach that separates Tx and Rx
% cirucuits, allow to identify the impact of each part separately: Tx,
% medium, Rx
% Here we provide a model using complete Sommerfeld-type expressions that
% can be accurately evaluated using a desktop computer
%
% As an example of the potential use of the present model we consider WUSN
% application:
%
% *(Vuran. Internet of Underground things in Precision Agriculture: Architecture and Technology Aspects). 
% -IOUT pplications are usually buried insum-meter layer. Our contribution: 0.4 m and 3 meters for other applications.
%  Buried Sensors are buried 
% -Dry soil Er 2 and 6 and Sigma 5-15 10^-4 and 10^-5 S/m
% *Simulation: consider upper soil of:
% *If we consider antennas can be arbitrarilly oriented: vertical, horizontal, 45 deg
% *Wet soil 70 mS/m er 12
% *Range: 100 dB of path loss
% * Sensors' location: UG2AG and UG2UG
% frequency sweep: 1kHz to 10MHz
% central layer 1m
%
% 
%
% Creation date: 1st February 2019
% Last modifiation:
%
% Author: PhD. Natalia Ayuso
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear all
clc

computations=0;
saving=1;

dir_path= '../../MATLAB_FEM_DATA/';
header=[dir_path 'UWSN_agriculture'];

filemat=[header '.mat'];

if computations==1   
global mu0 e0 c

e0=8.85e-12;
mu0=4*pi*1e-7;
c=1/sqrt(e0*mu0);

saving=1;           %1 for saving computations and 0 for avoit saving
% Emitter
he=[-3 -0.4 -0.1];
m0=1;

%f=logspace(3,7,300);     %frequency range 
f=logspace(3,7,10);     %frequency range 

% Media parameters
sigma1=70e-3;
er1=12;
e1=er1*e0;

sigma2=5e-4;
er2=2;
e2=er2*e0;

d=1;        %central layer width

% Receiver
hr=he;
rho=0;

py=rho;
px=0;

AbsTol=1e-16;
RelTol=0;

[Mf,Mhe]=meshgrid(f,he);

%AbsTol Default value 1e-10 (double)
%RelTol Default value 1e-6 (double)
%quadgk attempts to satisfy errbnd <= max(AbsTol,RelTol*|Q|). 
%This is absolute error control when |Q| is sufficiently small and relative error control when |Q| is larger. 
%For pure absolute error control use 'AbsTol' > 0 and'RelTol'= 0. 
%For pure relative error control use 'AbsTol' = 0. 
%Except when using pure absolute error control, the minimum relative tolerance is 'RelTol' >= 100*eps(class(Q)).

%%%%%%%
% HMD %
%%%%%%%

% fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nCALCULUS\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
normalization=0; 
factor_norm=1;
ULIM=Inf;

AbsTolx=AbsTol;
AbsToly=AbsTol;
AbsTolz=AbsTol;        
c_ref='cartesian';
aprox=0;
normalization=0;

toc_iter_HMD=[];
tic;
[index_he,index_f]=size(Mhe);

for i=1:index_he %he swept
    fprintf('f swept for he=%2.1f \n',he(i));
        tic;
        for j=1:index_f %f swept
            HTOL=integral_HMD_3_layers(m0,Mf(i,j),Mhe(i,j),0,5,Mhe(i,j),d,sigma1,sigma2,e1,e2,aprox,ULIM,ULIM,ULIM,AbsTolx,AbsToly,AbsTolz,RelTol,RelTol,RelTol,c_ref,normalization);
            AbsTolx=min(abs(HTOL(1))/10,AbsTol);
            AbsToly=min(abs(HTOL(2))/10,AbsTol);
            AbsTolz=min(abs(HTOL(3))/10,AbsTol);
            Hfree=VMDilhi (m0,Mf(i,j),Mhe(i,j),Mhe(i,j)+5,0,sigma1,e1,0,c_ref);
            Hfree0=VMDilhi (m0,Mf(i,j),Mhe(i,j),Mhe(i,j)+5,0,0,e0,0,c_ref);
            H_x(i,j) = HTOL(1); 
            H_y(i,j) = HTOL(2); 
            H_z(i,j) = HTOL(3); 
            H_x_free(i,j) = Hfree(1); 
            H_y_free(i,j) = Hfree(2); 
            H_z_free(i,j) = Hfree(3);
            H_x_free0(i,j) = Hfree0(1); 
            H_y_free0(i,j) = Hfree0(2); 
            H_z_free0(i,j) = Hfree0(3);
            j
        end
        H_HMD(i,:)=sqrt(abs(H_x(i,:)).^2+abs(H_y(i,:)).^2+abs(H_z(i,:)).^2);
        Hfree_total(i,:)=sqrt(abs(H_x_free(i,:)).^2+abs(H_y_free(i,:)).^2+abs(H_z_free(i,:)).^2);
        Hfree_total0(i,:)=sqrt(abs(H_x_free0(i,:)).^2+abs(H_y_free0(i,:)).^2+abs(H_z_free0(i,:)).^2);


        toc_iter_HMD=[toc_iter_HMD toc];  
end
    if saving
        save(filemat,'f','he','hr','Mf','Mhe','sigma1','er1','sigma2','er2','d','H_HMD','Hfree_total0','Hfree_total','toc_iter_HMD');
    end
else
    load (filemat)
end

figure(1)
loglog(f,(Hfree_total0(1,:)/Hfree_total0(1,1)),'LineWidth',1,'DisplayName','Free-space')
hold on
loglog(f,(Hfree_total(1,:)/Hfree_total0(1,1)),'LineWidth',1,'DisplayName','Infinite-Medium')
loglog(f,(H_HMD(3,:)/Hfree_total0(1,1)),'LineWidth',1,'DisplayName',['he=hr= ', num2str(he(3)),' m'])
loglog(f,(H_HMD(2,:)/Hfree_total0(1,1)),'LineWidth',1,'DisplayName',['he=hr= ', num2str(he(2)),' m'])
loglog(f,(H_HMD(1,:)/Hfree_total0(1,1)),'LineWidth',1,'DisplayName',['he=hr= ', num2str(he(1)),' m'])
xlabel('f (Hz)')
ylabel('|H_{NORMALIZED}| (dB). Dipole separation = 5 m')
legend('show')

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nRESULTS CALCULUS\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')

fprintf('Total time VMD %s hours %s min \n', num2str(round(sum(toc_iter_HMD/3600))),num2str(round(sum(toc_iter_HMD/60))));
