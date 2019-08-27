%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAGNETIC FIELD SOMMERFELD INTEGRALS DUE TO A VMD
% 3-LAYERS MODEL: AIR-SOIL1-SOIL2
%
%   z
%   ^  
%   |
%   | ·(px, py, hr)        AIR   
%   -----------------------------------> rho
%   |                               ^
%   |                               |
%   |                               |
%   | * (he)        MEDIUM 1    d
%   |                               |
%   |                               |
%   |                               v
%   ----------------------------------
%   |
%   |                   MEDIUM 2
%
%
% Observation point: px, py, py
% Source: he, hr
% Axial symmetry: cylindrical coordinate system
%
%
% AUTHOR: PhD. Natalia Ayuso Escuer at nayuso@unizar.es
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES DESCRIPTION
%
% m0: magnetic moment
% f: frequency
% he: z coordiante-source
% hr: z coordinate-observation point
% offset: rho coordinate-observation point
% d: width layer 1
% sigma1: medium 1 conductivity
% sigma2: medium 2 conductivity
% e1: medium 1 relative permittivity
% e2: medium 2 relative permittivity
% disp_current: 1--> displacement currents neglected; 
% ULIMr: upper integration limit for z (or rho) component
% TOLABSr: Absolute Tolerance for r (or z) component
% TOLRELr: Relative Tolerance for r (or z) component
% c_ref: Coordinate system for magnetic field expressions
% normalization: 1-->change of integration variable, 0-->h=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H3num = VMD_3layers (m0,f,he,hr,px,py,d,sigma1,sigma2,e1,e2,disp_current,ULIMrho, ULIMz, TOLABSr,TOLABSz,TOLRELr,TOLRELz,c_ref,normalization)

global Hr_c3 Hz_c3 valhr valrho mu0 e0

offset = sqrt(px.^2+py.^2);
rphi = atan2(py,px);                %used in cylindric coordinate system
rr = offset;
rz = hr-he;

r = sqrt(px.^2+py.^2+rz.^2);        %r: tx-rx distance


%Used for normalization
if he<0
    h=abs(he);
elseif hr<0
    h = abs(hr);
else
    h = 1;
end

if h<1
    h=1;
end

if normalization==0
    h=1;
end

%PRIMARY FIELDS
j = sqrt(-1);
w = 2*pi*f;
ec1 = e1 - j*(sigma1/w); %complex permittivity of med1
ec2 = e2 - j*(sigma2/w); %complex permittivity of med1


if disp_current == 1  %if displacement currents are neglected
    kc0 = 0;
    kc1 = sqrt(-j*mu0*sigma1*w);
    kc2 = sqrt(-j*mu0*sigma2*w);
else
    kc0 = w*sqrt(mu0*e0);
    kc1 = w*sqrt(mu0*ec1);
    kc2 = w*sqrt(mu0*ec2);
end

if hr >= 0
    kc = kc0;
elseif hr < (-d)
    kc = kc2;
else
    kc = kc1;
end

R = sqrt(px.^2+py.^2+(hr-he).^2);

%PRIMARY FIELDS
Hr_c1V = - (m0/(4*pi)).*offset.*(hr-he).*((kc.^2./(R.^2))-(3*j*kc./(R.^3))-(3./(R.^4))).*(exp(-j*kc.*R))./R;
Hphi_c1V = zeros(size(Hr_c1V));
Hz_c1V =  (m0/(4*pi)).*((kc.^2 - (j*kc./R) - (1./(R.^2))) - ((hr-he).^2).*((kc.^2./(R.^2))-(3*j*kc./(R.^3))-(3./(R.^4)))).*(exp(-j*kc.*R))./R;

%SECONDARY FIELDS
Hr_c3 = (m0/(4*pi*(abs(h)^3)))*integral (@(x)VMD_integrand_3layers(x,f,he,hr,offset,d,sigma1,sigma2,e1,e2,a,0,disp_current,normalization),0,ULIMr,'AbsTol',TOLABSr,'Reltol',TOLRELr);

Hphi_c3 = zeros(size(Hr_c3));

Hz_c3 = (m0/(4*pi*(abs(h)^3)))*integral (@(x)VMD_integrand_3layers(x,f,he,hr,offset,d,sigma1,sigma2,e1,e2,a,1,disp_current,normalization),0,ULIMz,'AbsTol',TOLABSz,'Reltol',TOLRELz);


%SOURCE AND OBSERVER IN THE SAME MEDIUM: PRIMARY + SECONDARY
%Tx-Rx in air=(he >= 0)&(hr >=0))
%Tx-Rx in medium 2=((he <= -d)&(hr <= -d))
%Tx-Rx in medium 1=(hr <= -d))|(((-d < he)&(he < 0))&((-d < hr)&(hr < 0))
if ((he >= 0)&(hr >=0))|((he < -d)&(hr < -d))|(((-d <= he)&(he < 0))&((-d <= hr)&(hr < 0)))
    Hr_c3  = Hr_c3 + Hr_c1V;
    Hphi_c3= Hphi_c3 + Hphi_c1V;
    Hz_c3  = Hz_c3 + Hz_c1V;
end

switch c_ref
    case 'cartesian'
        phi = rphi;
        Hx_k3 = cos(phi).*Hr_c3 - sin(phi).*Hphi_c3;
        Hy_k3 = sin(phi).*Hr_c3 + cos(phi).*Hphi_c3;
        Hz_k3 = Hz_c3;
        H3num = [Hx_k3;Hy_k3;Hz_k3]; 
    case 'cylindric'
        H3num = [Hr_c3;Hphi_c3;Hz_c3]; 
    case 'spheric'
        if ((hr == 0)&(px ==0)& (py == 0))
            tita = 0;
        else
            tita = acos(hr./sqrt(px.^2+py.^2+hr.^2));
        end
        Hr_e3 = sin(tita).*Hr_c3 + cos(tita).*Hz_c3;
        Htita_e3 = cos(tita).*Hr_c3 - sin(tita).*Hz_c3;
        Hphi_e3 = Hphi_c3;
        H3num= [Hr_e3; Htita_e3; Hz_c3];
end

function Hint=VMD_integrand_3layers(x,f,he,hr,offset,d,sigma1,sigma2,e1,e2,a,componente,disp_current,normalization

%%%%%%%%%%%%%
% CONSTANTS %
%%%%%%%%%%%%%
global mu0 e0;

if he<0
    h=abs(he);
elseif hr<0
    h = abs(hr);
else
    h = 1;
end

%to avoid problems with convergence
if h<1
    h=1;
end

if normalization==0
    h=1;
end

j = sqrt(-1);
w = 2*pi*f;
z = hr;         %z-coordiante of receiver 
rho = offset;   %rho-coordinate of recevier

if disp_current == 1
    gamma0 = 0;
    gamma1 = sqrt(j*sigma1*mu0*w); 
    gamma2 = sqrt(j*sigma2*mu0*w); 
else
    gamma0 = sqrt(-(w^2*e0*mu0)); 
    gamma1 = sqrt((j*sigma1*mu0*w)-(w^2*e1*mu0)); 
    gamma2 = sqrt((j*sigma2*mu0*w)-(w^2*e2*mu0)); 
end


%TX-RX's mediums
emitter = 0; %EMITTER MEDIUM: 0 air of air-med1 interface; 1 med 1 or med1-2 interface; 2 med 2
if ((he < 0)&(he >= -d)) 
    emitter = 1;
elseif (he < -d)
    emitter = 2;
end        
medium = 0;  %RECEIVER MEDIUM: 0 air of air-med1 interface; 1 med 1; 2 med 2
if ((z < 0)&(z >= -d)) 
    medium = 1;
elseif z < -d
    medium = 2;
end

radial = 0; %radial=1: for radial component and radial=0 for axial component

%%%%%%%%%%%%%
% INTEGRAND %
%%%%%%%%%%%%%

%BESSEL FUNCTIONS
J0 = besselj(0,(x/h)*rho); 
J1 = besselj(1,(x/h)*rho); 

u0 = sqrt((x/h).^2 + gamma0.^2); %propagation constant air,
u1 = sqrt((x/h).^2 + gamma1.^2); %medium 1,
u2 = sqrt((x/h).^2 + gamma2.^2); %medium 2

%%%%%%%%%
% T, R  %
%%%%%%%%%
R01h = (u0-u1)./(u0+u1);
R02h = (u0-u2)./(u0+u2);
R10h = (u1-u0)./(u1+u0);
R12h = (u1-u2)./(u1+u2);
R20h = (u2-u0)./(u2+u0);
R21h = (u2-u1)./(u2+u1);
T01h = 2*u0./(u0+u1);
T21h = 2*u2./(u1+u2);


if emitter == 1 %Source med 1
    C1 = ((x/h)./u1).*R12h.*(exp(-u1.*(z+2*d+he)) + R10h.*exp(-u1.*(z+2*d-he)))./(1 - R10h.*R12h.*exp(-2*u1*d));
    B1 = ((x/h)./u1).*R10h.*(exp(u1.*(z+he)) + R12h.*exp(u1.*(z-2*d-he)))./(1 - R10h.*R12h.*exp(-2*u1*d));    
    A1 = (((x/h)./u1).*exp(-u1*abs(he)-u0.*z))+ ((x/h)./u1).*R10h.*(exp(u1*he-u0.*z) + R12h.*exp(-u1*(2*d+he)-u0.*z))./(1 - R10h.*R12h.*exp(-2*u1*d)) + ((x/h)./u1).*R12h.*(exp(-u1*(2*d+he)-u0.*z) + R10h.*exp(-u1*(2*d-he)-u0.*z))./(1 - R10h.*R12h.*exp(-2*u1*d));
    D1 = (((x/h)./u1).*exp(-u1*(d-abs(he))+u2*(d+z)) + ((x/h)./u1).*R10h.*(exp(-u1*(d-he)+u2*(d+z)) + R12h.*exp(-u1*(3*d+he)+u2*(d+z)))./(1 - R10h.*R12h.*exp(-2*u1*d)) + ((x/h)./u1).*R12h.*(exp(-u1*(d+he)+u2*(d+z)) + R10h.*exp(u1*(he-d)+u2*(d+z)))./(1 - R10h.*R12h.*exp(-2*u1*d)));    

    if medium == 0 %Air or air-med1 interface

        if componente == 0;
            Hint = A1.*(x*h).*u0.*J1; %Radial
        else
            Hint = A1.*(x.^2).*J0; %Axial 
        end

    elseif medium == 1 %Med 1 or med1-med2 interface
        if componente == 0;
              Hint = x.*h.*u1.*(-B1+C1).*J1;    %Radial 
        else
            Hint = (x.^2).*(B1+C1).*J0;         %Axial  
        end
    else %Med2
        if componente == 0;
            Hint = -D1.*(x*h).*u2.*J1;          %Radial 
        else
            Hint = D1.*(x.^2).*J0;              %Axial
        end
    end

elseif emitter == 0 %Source Air
    
        A0=-((x/h)./u0).*exp(-u0*(he+z))+(1+R10h).*T01h.*R12h.*((x/h)./u0).*exp(-2*u1*d-u0*(z+he))./(1 - R10h.*R12h.*exp(-2*u1*d))+T01h.*((x/h)./u0).*exp(-u0*(z+he))./(1 - R10h.*R12h.*exp(-2*u1*d));
        C0 = T01h.*R12h.*((x/h)./u0).*exp(u1*(-2*d+z)-u0*he)./(1 - R10h.*R12h.*exp(-2*u1*d));
        D0=(R10h.*exp(-2*u1*d)+1).*(R12h.*T01h.*((x/h)./u0).*exp(-u1*d-u0*he+u2*(d+z))./(1 - R10h.*R12h.*exp(-2*u1*d)))+T01h.*((x/h)./u0).*exp(-u0*he-u1*d+u2*(d+z));
    
    if medium == 0 %Air or air-med1 interface
        if componente == 0;        
            Hint = x*h.*u0.*A0.*J1; %Radial
       else
            Hint = (x.^2).*A0.*J0;  %Axial 
        end
    elseif medium == 1 %Med1 or med1-med2 interface
        if componente == 0; 
           Hint = x.*h.*u1.*(-C0.*R10h+C02-T01h.*((x/h)./u0).*exp(u1*z-u0*he)).*J1; %Radial 

        else
            Hint = (x.^2).*(C0.*(R10h+1)+T01h.*((x/h)./u0).*exp(u1*z-u0*he)).*J1;   %Axial

        end
    else %Med2
        if componente == 0;
            Hint = D0.*(-x*h).*u2.*J1;  %Radial        
        else
            Hint = D0.*(x.^2).*J0;      %Axial
        end
    end

elseif emitter == 2%Source med2
   
    C2 = (((x/h)./u2).*exp(-u2*abs(d-abs(he))).*exp(-u1.*(d+z)).*T21h)./(1-R10h.*R12h.*exp(-2*u1*d));
    A2 = C2.*exp((u1-u0).*z).*(R10h + 1);
    D2 = ((x/h)./u2).*exp(u2*(2*d+z-abs(he))).*(-1+(T21h.*(R10h.*exp(-2*u1*d)+1)./(1-R10h.*R12h.*exp(-2*u1*d))));
 
    if medium == 0 %Air or air-med1 interface
        if componente == 0; 
            Hint = A2.*(x*h).*u0.*J1;   %Radial
        else
            Hint = A2.*(x.^2).*J0;      %Axial 
        end   
  elseif medium == 1 %Med 1 or med1-med2 interface
        if componente == 0; 
            Hint = x.*h.*u1.*C2.*(-R10h.*exp(u1.*z)+1).*J1; %Radial 
        else
            Hint = (x.^2).*C2.*(R10h.*exp(u1.*z)+1).*J0;    %Axial 
        end
  else %Med2
        if componente == 0;        
            Hint = -x*h.*u2.*D2.*J1; %Radial
        else
            Hint = (x.^2).*D2.*J0;  %Axial
        end
  end
end