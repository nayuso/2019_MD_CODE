%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAGNETIC FIELD INFINITE MEDIUM
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
% sigma: conductivity
% e: relative permittivity
% disp_current: 1--> displacement currents neglected; 
% c_ref: Coordinate system for magnetic field expressions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function H1V = VMDilhi (m0,f,he,hr,offset,sigma,e,disp_current,c_ref)

%%%%%%%%%%%%
%CONSTANTS %
%%%%%%%%%%%%
global mu0 e0;

j = sqrt(-1);
w = 2*pi*f;
ec = e - j*(sigma/w); %permitividad compleja

if disp_current == 1
    kc = sqrt(-j*mu0*sigma*w); 
else
    kc = w*sqrt(mu0*ec); 
end

R = sqrt(offset.^2+(hr-he).^2); 
phi = 0; 

if ((hr == 0)&&(offset ==0))
    tita = 0;
else
    tita = acos(hr./sqrt(offset.^2+hr.^2));
end



%Cylindric
Hr_c1V = - (m0/(4*pi)).*offset.*(hr-he).*((kc.^2./(R.^2))-(3*j*kc./(R.^3))-(3./(R.^4))).*(exp(-j*kc.*R))./R;
Hphi_c1V = zeros(size(Hr_c1V));
Hz_c1V =  (m0/(4*pi)).*((kc.^2 - (j*kc./R) - (1./(R.^2))) - ((hr-he).^2).*((kc.^2./(R.^2))-(3*j*kc./(R.^3))-(3./(R.^4)))).*(exp(-j*kc.*R))./R;


switch c_ref
    case 'cartesian'
        Hx_k1V = cos(phi).*Hr_c1V - sin(phi).*Hphi_c1V;
        Hy_k1V = sin(phi).*Hr_c1V + cos(phi).*Hphi_c1V;
        Hz_k1V = Hz_c1V;
        H1V = [Hx_k1V; Hy_k1V; Hz_k1V];     %coordenadas cartesianas
    case 'cylindric'
        H1V = [Hr_c1V; Hphi_c1V; Hz_c1V];   %coordenadas cylindric
    case 'spheric'
        Hr_e1V = sin(tita).*Hr_c1V + cos(tita).*Hz_c1V;
        Htita_e1V = cos(tita).*Hr_c1V - sin(tita).*Hz_c1V;
        Hphi_e1V = Hphi_c1V;
        H1V = [Hr_e1V; Htita_e1V; Hphi_e1V];     %spheric coordinates
end
