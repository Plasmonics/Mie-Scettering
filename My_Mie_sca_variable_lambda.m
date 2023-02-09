% Mie scattering
%% ------------------------------------
% Mie scattering cross-section for wavelength dependent refractive index
%  Last update: 10/14/2018
%  Author:Hosna Sultana

% Acknowledgment and reference: http://www.guillaume.baffou.com/pdf/Mie_Theory.pdf
%% ------------------------------------

% display_by_Mie = 1; % 500 nm lambda for Water droplet of a = 272.58467103348891669 nm
%   display_by_Mie = 2; % 200 -1000 nm lambda for Water droplet of a = 272.58467103348891669 nm
%  display_by_Mie = 3; % 200 -1000 nm lambda for Ag NP of a = 5 nm
   display_by_Mie = 4; % 200 -1000 nm lambda n = n(lambda)for Ag NP of a = 5 nm Winsemius
%  display_by_Mie = 5; % 200 -1000 nm lambda n = n(lambda)for Ag NP of a = 5 nm, Babar 
% display_by_Mie = 6; % 200 -1000 nm lambda for Ag NP of a = 5 nm
% display_by_Mie = 7; % 200 -1000 nm lambda for Ag NP of a = 5 nm



%% 0.5 micrometer lambda for Water droplet of a = 0.2725846710334890 micrometer
if display_by_Mie == 1
    
m_p = (1.3660 + 0.005i);
% m_p = (0.135 + 3.99i);
lambda = 0.550;
k = (2* pi /lambda);
% r_p = 50; %in nm
% rho = k * r_p;
a = 0.272584671033489;   % nm
x = k * a;
z = m_p * x;
n_max = round(x + 4.05 *(x) .^(1/3) + 2);
m_m = 1;
n = (1:n_max);
%spherical bessel (n,z) = besselj(n+1/2,z)*sqrt(pi/(2*z))
%shi_n(x) = rho * j_n(x)

prefacx = x .*(pi./ (2*x)).^(1/2);
prefacz = z .*(pi./ (2*z)).^(1/2);
shi_x = prefacx .* besselj(n+0.5,x);
kai_x = prefacx .* (besselj(n+0.5,x)+1i*bessely(n+ 0.5,x));
shi_z = prefacz .* besselj(n+0.5,z);
y0 = prefacx .* bessely (n + 0.5, x);
y1 = [-cos(x), y0(1 : n_max-1)];
% d(shi_x) = shi_x(n-1,x)- n*shi_x(n,x)/x
% let,   shi_x(n-1,x) = shi_xn
% let,   shi_z(n-1,z) = shi_zn

shi_xn = [sin(x), shi_x(1 : n_max-1)];
shi_zn = [sin(z), shi_z(1 : n_max-1)];
dshi_x = [(shi_xn - n/x.*shi_x)] ;
%  d(kai) =  kai(n-1,x)- n*kai(n,x)/x. ;
dkai_x =  (shi_xn + 1i * y1)- n/x.*(shi_x + 1i*y0);
dshi_z =  (shi_zn - n/z.*shi_z);

a_n = (m_p * shi_z .* dshi_x - shi_x .* dshi_z)./(m_p * shi_z .* dkai_x - kai_x .* dshi_z);
b_n = (shi_z .* dshi_x - m_p * shi_x .* dshi_z)./(shi_z .* dkai_x - m_p * kai_x .* dshi_z);

sigma_sca = 2*pi./ k .^2 .* sum ((2*n+1) .* (abs(a_n) .* abs(a_n)  +  abs(b_n) .* abs(b_n)));
sigma_ext = 2*pi./ k .^2 .* sum ((2*n+1) .* real(a_n + b_n));
sigma_abso = (sigma_ext - sigma_sca);
Q_sca = sigma_sca/(pi * a^2);
Q_ext = sigma_ext/(pi * a^2);
Q_abso = sigma_abso/(pi * a^2);

ZP  = [shi_z length(n_max)./shi_z length(n_max -1)];
ZA = dshi_z ./shi_z ; 





elseif display_by_Mie == 2
%% 0.2 - 1 micrometer lambda for Water droplet of a = 0.2725846710334890 micrometer
  m_p = (1.3660 + 0.005i);  
%     m_p = (0.135 + 3.99i);
m_m = 1;
%   lambda = 0.550;
   lambda =200: 1: 1000;
   L= 200: 1: 1000;% in nm
i=1;
for lambda =200: 1: 1000
k = (2* pi ./lambda);
% a = 5;   % nm
a = 0.2725846710334890*10^3; 
x = (a .* k);
z = (m_p .* x);
 n_max = round(x + 4.05 .*(x) .^(1/3) + 2);
 n = 1:n_max; 
%spherical bessel (n,z) = besselj(n+1/2,z)*sqrt(pi/(2*z))
%shi_n(x) = rho * j_n(x)
prefacx = x .*(pi./ (2*x)).^(1/2);
prefacz = z .*(pi./ (2*z)).^(1/2);
shi_x = prefacx .* besselj(n+0.5,x);
kai_x = prefacx .* (besselj(n+0.5,x)+1i*bessely(n+ 0.5,x));
shi_z = prefacz .* besselj(n+0.5,z);
y0 = prefacx .* bessely (n + 0.5, x);
y1 = [-cos(x), y0(1 : n_max-1)];

shi_xn = [sin(x), shi_x(1 : n_max-1)];
shi_zn = [sin(z), shi_z(1 : n_max-1)];

dshi_x = [(shi_xn - n/x.*shi_x)] ;
%  d(kai) =  kai(n-1,x)- n*kai(n,x)/x. ;
dkai_x =  (shi_xn + 1i * y1)- n./x.*(shi_x + 1i*y0);
dshi_z =  (shi_zn - n/z.*shi_z);

a_n = (m_p .* shi_z .* dshi_x - shi_x .* dshi_z)./(m_p .* shi_z .* dkai_x - kai_x .* dshi_z);
b_n = (shi_z .* dshi_x - m_p .* shi_x .* dshi_z)./(shi_z .* dkai_x - m_p .* kai_x .* dshi_z);
 
sigma_sca(i) = 2*pi./ k .^2 .* (sum ((2*n+1) .* (abs(a_n) .* abs(a_n)  +  abs(b_n) .* abs(b_n))));
sigma_ext(i) = 2*pi./ k .^2 .* sum ((2*n+1) .* real(a_n + b_n));
sigma_abso(i) = (sigma_ext(i) - sigma_sca(i));
Q_sca(i) = sigma_sca(i)/(pi * a^2);
Q_ext(i) = sigma_ext(i)/(pi * a^2);
Q_abso(i) = sigma_abso(i)/(pi * a^2);
lambda = lambda + 1;
i = i+1;
end

ZP  = [shi_z length(n_max)./shi_z length(n_max -1)];
ZA = dshi_z ./shi_z ;  

fontsize = 14;
linewidth = 1.5;
 markersize = 1;
 
  subplot(2,1,1)
plot(L, Q_sca,'.-','color',[1.0 0.564 0.141],'LineWidth',linewidth,'MarkerSize',markersize);
  hold on
plot(L, Q_ext,'c-','LineWidth',linewidth,'MarkerSize',markersize);
plot(L, Q_abso,'.-','color',[0.0 0.4 0.0],'LineWidth',linewidth,'MarkerSize',markersize);
  hold off 
  legend('Q_{sca}','Q_{ext}','Q_{abso}','All','Location','northeast');
grid on
    set(gca,'FontSize',fontsize)
    set(gcf,'Color','w')
    ylabel('Q_{sca},Q_{ext},Q_{abso} ')
    xlabel('\lambda (nm)')  
    
 subplot(2,1,2)
plot(L, sigma_sca,'.-','color',[0.7 0.8 0.3],'LineWidth',linewidth,'MarkerSize',markersize);
  hold on
plot(L, sigma_ext,'.-','color',[0.6 0.1 1.0],'LineWidth',linewidth,'MarkerSize',markersize);
plot(L, sigma_abso,'.-','color',[0.0 0.4 0.5],'LineWidth',linewidth,'MarkerSize',markersize);
  hold off 
  legend('\sigma_{sca}','\sigma_{ext}','\sigma_{abso}','All','Location','northeast');
grid on
    set(gca,'FontSize',fontsize)
    set(gcf,'Color','w')
    ylabel('\sigma_{sca},\sigma_{ext},\sigma_{abso} in mm^2')
    xlabel('\lambda (nm)')     
 annotation('textbox',[0.34 0.939 0.7 0.024],'String',{'Mie scattering for Water droplet of radius 272.5845 nm of optical index (1.3660 + 0.005i)'},'FitBoxToText','on');   
% annotation('textbox',[0.34 0.939 0.7 0.024],'String',{'Mie scattering for Silver NP radius 272.5845 nm of optical index (0.135 + 3.99i)'},'FitBoxToText','on');   

%  
%  I = I_o exp(-(b_sca + b_abs)l)   for number of particles
% b_x = <sigm_x> (particle number/m^3)





elseif display_by_Mie == 3
%% 200 -1000 nm lambda for Ag NP of a = 5 nm
%  m_p = (1.3660 + 0.005i);  
   m_p = (0.135 + 3.99i);
m_m = 1;
%   lambda = 0.550;
   lambda =200: 1: 1000;
   L= 200: 1: 1000;% in nm
i=1;
for lambda =200: 1: 1000
k = (2* pi ./lambda);
 a = 5;   % nm
% a = 0.2725846710334890; 
x = (a .* k);
z = (m_p .* x);
 n_max = round(x + 4.05 .*(x) .^(1/3) + 2);
 n = 1:n_max; 
%spherical bessel (n,z) = besselj(n+1/2,z)*sqrt(pi/(2*z))
%shi_n(x) = rho * j_n(x)
prefacx = x .*(pi./ (2*x)).^(1/2);
prefacz = z .*(pi./ (2*z)).^(1/2);
shi_x = prefacx .* besselj(n+0.5,x);
kai_x = prefacx .* (besselj(n+0.5,x)+1i*bessely(n+ 0.5,x));
shi_z = prefacz .* besselj(n+0.5,z);
y0 = prefacx .* bessely (n + 0.5, x);
y1 = [-cos(x), y0(1 : n_max-1)];

shi_xn = [sin(x), shi_x(1 : n_max-1)];
shi_zn = [sin(z), shi_z(1 : n_max-1)];

dshi_x = [(shi_xn - n/x.*shi_x)] ;
%  d(kai) =  kai(n-1,x)- n*kai(n,x)/x. ;
dkai_x =  (shi_xn + 1i * y1)- n./x.*(shi_x + 1i*y0);
dshi_z =  (shi_zn - n/z.*shi_z);

a_n = (m_p .* shi_z .* dshi_x - shi_x .* dshi_z)./(m_p .* shi_z .* dkai_x - kai_x .* dshi_z);
b_n = (shi_z .* dshi_x - m_p .* shi_x .* dshi_z)./(shi_z .* dkai_x - m_p .* kai_x .* dshi_z);
 
sigma_sca(i) = 2*pi./ k .^2 .* (sum ((2*n+1) .* (abs(a_n) .* abs(a_n)  +  abs(b_n) .* abs(b_n))));
sigma_ext(i) = 2*pi./ k .^2 .* sum ((2*n+1) .* real(a_n + b_n));
sigma_abso(i) = (sigma_ext(i) - sigma_sca(i));
Q_sca(i) = sigma_sca(i)/(pi * a^2);
Q_ext(i) = sigma_ext(i)/(pi * a^2);
Q_abso(i) = sigma_abso(i)/(pi * a^2);
lambda = lambda + 1;
i = i+1;
end

ZP  = [shi_z length(n_max)./shi_z length(n_max -1)];
ZA = dshi_z ./shi_z ;  


 fontsize = 14;
linewidth = 1.5;
 markersize = 1;
 
  subplot(2,1,1)
plot(L, Q_sca,'.-','color',[1.0 0.564 0.141],'LineWidth',linewidth,'MarkerSize',markersize);
  hold on
plot(L, Q_ext,'c-','LineWidth',linewidth,'MarkerSize',markersize);
plot(L, Q_abso,'.-','color',[0.0 0.4 0.0],'LineWidth',linewidth,'MarkerSize',markersize);
  hold off 
  legend('Q_{sca}','Q_{ext}','Q_{abso}','All','Location','northeast');
grid on
    set(gca,'FontSize',fontsize)
    set(gcf,'Color','w')
    ylabel('Q_{sca},Q_{ext},Q_{abso} ')
    xlabel('\lambda (nm)')  
    
 subplot(2,1,2)
plot(L, sigma_sca,'.-','color',[0.7 0.8 0.3],'LineWidth',linewidth,'MarkerSize',markersize);
  hold on
plot(L, sigma_ext,'.-','color',[0.6 0.1 1.0],'LineWidth',linewidth,'MarkerSize',markersize);
plot(L, sigma_abso,'.-','color',[0.0 0.4 0.5],'LineWidth',linewidth,'MarkerSize',markersize);
  hold off 
  legend('\sigma_{sca}','\sigma_{ext}','\sigma_{abso}','All','Location','northeast');
grid on
    set(gca,'FontSize',fontsize)
    set(gcf,'Color','w')
    ylabel('\sigma_{sca},\sigma_{ext},\sigma_{abso} in mm^2')
    xlabel('\lambda (nm)')     
 annotation('textbox',[0.34 0.939 0.7 0.024],'String',{'Mie scattering for Silver NP radius  5 nm of optical index (0.135 + 3.99i)'},'FitBoxToText','on');   


 
%  I = I_o exp(-(b_sca + b_abs)l)   for number of particles
% b_x = <sigm_x> (particle number/m^3)



elseif display_by_Mie == 4
%% 200 -1000 nm lambda for Ag NP of a = 5 nm with variable optical index

op_index  = xlsread('Ag_n_k_P_Winsemius_Kampen_1976.xlsx');
energy_eV_for_lambda= op_index(:,1);

Wavelength = op_index(:,1);
Ag_n_= op_index(:,2);
Ag_kappa_= op_index(:,3);
% Ag_reflectance= op_index(:,5);

interploation_Points = 10^6;

wavelength = zeros(interploation_Points, 1);
Ag_n = zeros(interploation_Points, 1);
Ag_kappa = zeros(interploation_Points, 1);
% Ag_R = zeros(interploation_Points, 1);

    Wavelength(122:end, :) = [];
    wv = Wavelength;
    Agn = op_index(:,2);
    Agkappa = op_index(:,3);
%     AgR = op_index(:,5);
  wv(isinf(wv)) = nan;
  fillmissing(wv, 'linear'); 
  wv = fillmissing(wv, 'linear');


   X = linspace(min(wv),max(wv),interploation_Points);
   [C,~,idwv] = unique(wv,'stable');
   val1 = accumarray(idwv,Agn,[],@mean);
   val2 = accumarray(idwv,Agkappa,[],@mean);
%    val3 = accumarray(idwv,AgR,[],@mean);
   Y = interp1(C,val1,X,'linear','extrap'); 
   Z = interp1(C,val2,X,'linear','extrap');
%    V = interp1(C,val3,X,'linear','extrap');
% wavelength = wavelength + X';
% Ag_n = Ag_n + Y';
wavelength = X;
Ag_n = Y;
Ag_kappa = Z;  
% Ag_R = V;


M_p = (Ag_n + (Ag_kappa).*1i);
m_p = linspace(min(M_p),max(M_p),interploation_Points);
%   m_p = (0.135 + 3.99i);
m_m = 1;
   Lambda = linspace(min(wavelength),max(wavelength),interploation_Points);
%    lambda =200: 1: 1000;
   L= Lambda;% in nm
for KK = 1 : length(L)
  this_L = L(KK);
  this_m_p = m_p(KK);

k = (2* pi ./this_L);
 a = 5;   % nm
% a = 0.2725846710334890; 
x = (a .* k);
z = (a.*this_m_p .* 2.*pi./this_L);

  n_max = round(x + 4.05 .*(x) .^(1/3) + 2);
 n = 1:n_max; 
%spherical bessel (n,z) = besselj(n+1/2,z)*sqrt(pi/(2*z))
%shi_n(x) = rho * j_n(x)
prefacx = x .*(pi./ (2*x)).^(1/2);
prefacz = z .*(pi./ (2*z)).^(1/2);
shi_x = prefacx .* besselj(n+0.5,x);
kai_x = prefacx .* (besselj(n+0.5,x)+1i*bessely(n+ 0.5,x));
shi_z = prefacz .* besselj(n+0.5,z);
y0 = prefacx .* bessely (n + 0.5, x);
y1 = [-cos(x), y0(1 : n_max-1)];

shi_xn = [sin(x), shi_x(1 : n_max-1)];
shi_zn = [sin(z), shi_z(1 : n_max-1)];

dshi_x = [(shi_xn - n/x.*shi_x)] ;
%  d(kai) =  kai(n-1,x)- n*kai(n,x)/x. ;
dkai_x =  (shi_xn + 1i * y1)- n./x.*(shi_x + 1i*y0);
dshi_z =  (shi_zn - n/z.*shi_z);

 a_n = (this_m_p .* shi_z .* dshi_x - shi_x .* dshi_z)./(this_m_p .* shi_z .* dkai_x - kai_x .* dshi_z);
 b_n = (shi_z .* dshi_x - this_m_p .* shi_x .* dshi_z)./(shi_z .* dkai_x - this_m_p .* kai_x .* dshi_z);

 
sigma_sca(KK) = 2*pi./ k .^2 .* (sum ((2*n+1) .* (abs(a_n) .* abs(a_n)  +  abs(b_n) .* abs(b_n))));
sigma_ext(KK) = 2*pi./ k .^2 .* sum ((2*n+1) .* real(a_n + b_n));
sigma_abso(KK) = (sigma_ext(KK) - sigma_sca(KK));
Q_sca(KK) = sigma_sca(KK)/(pi * a^2);
Q_ext(KK) = sigma_ext(KK)/(pi * a^2);
Q_abso(KK) = sigma_abso(KK)/(pi * a^2);
 
this_L = this_L + 1;
this_m_p=this_m_p+1;
KK = KK+1;

ZP  = [shi_z length(n_max)./shi_z length(n_max -1)];
ZA = dshi_z ./shi_z ;  
end
 
fontsize = 14;
linewidth = 1.5;
 markersize = 1;
 
  subplot(3,1,1)
plot(L, Q_sca,'.-','color',[1.0 0.564 0.141],'LineWidth',linewidth,'MarkerSize',markersize);
  hold on
plot(L, Q_ext,'c-','LineWidth',linewidth,'MarkerSize',markersize);
plot(L, Q_abso,'.-','color',[0.0 0.4 0.0],'LineWidth',linewidth,'MarkerSize',markersize);
  hold off 
  legend('Q_{sca}','Q_{ext}','Q_{abso}','All','Location','northeast');
grid on
    set(gca,'FontSize',fontsize)
    set(gcf,'Color','w')
    ylabel('Q_{sca},Q_{ext},Q_{abso} ')
%     xlabel('\lambda (nm)')  
    xlim([200 1000])
    
     
 
 subplot(3,1,2)
plot(L, sigma_sca,'.-','color',[0.7 0.8 0.3],'LineWidth',linewidth,'MarkerSize',markersize);
  hold on
plot(L, sigma_ext,'.-','color',[0.6 0.1 1.0],'LineWidth',linewidth,'MarkerSize',markersize);
plot(L, sigma_abso,'.-','color',[0.0 0.4 0.5],'LineWidth',linewidth,'MarkerSize',markersize);
  hold off 
  legend('\sigma_{sca}','\sigma_{ext}','\sigma_{abso}','All','Location','northeast');
grid on
    set(gca,'FontSize',fontsize)
    set(gcf,'Color','w')
    ylabel('\sigma_{sca},\sigma_{ext},\sigma_{abso} in mm^2')
%     xlabel('\lambda (nm)')     
 annotation('textbox',[0.34 0.939 0.7 0.024],'String',{'Mie scattering for Silver NP radius  5 nm with optical index as a function of \lambda '},'FitBoxToText','on');   
  xlim([200 1000])
  
  subplot(3, 1, 3)
  plot(L, Ag_n,'g-','LineWidth',linewidth,'MarkerSize',markersize);
  hold on
  plot(L, Ag_kappa,'b-','LineWidth',linewidth,'MarkerSize',markersize);
%   plot(L, Ag_R,'m-','LineWidth',linewidth,'MarkerSize',markersize);
  hold off
  legend('n for Ag','\kappa for Ag','Reflectance','All','Location','southeast');
  grid on
   xlim([200 1000])
   ylabel('optical index of Ag')
   xlabel('\lambda (nm)')
   
   
  handaxes1 = axes('position',[0.62,0.77,0.12,0.12]);
      fontsize = 7;
     linewidth = 0.5;
%  markersize = 0.5;
plot(L, Q_sca,'.-','color',[1.0 0.564 0.141],'LineWidth',linewidth,'MarkerSize',markersize);
  hold on
plot(L, Q_ext,'c-','LineWidth',linewidth,'MarkerSize',markersize);
plot(L, Q_abso,'.-','color',[0.0 0.4 0.0],'LineWidth',linewidth,'MarkerSize',markersize);
  hold off 
  set(handaxes1, 'Box', 'off')
 grid on
 xlim([260 300])
  ylim([0.78 0.82])
 
  handaxes2 = axes('position',[0.62,0.48,0.12,0.12]);
      fontsize = 7;
     linewidth = 0.5;
%  markersize = 0.5;
plot(L, sigma_sca,'.-','color',[0.7 0.8 0.3],'LineWidth',linewidth,'MarkerSize',markersize);
  hold on
plot(L, sigma_ext,'.-','color',[0.6 0.1 1.0],'LineWidth',linewidth,'MarkerSize',markersize);
plot(L, sigma_abso,'.-','color',[0.0 0.4 0.5],'LineWidth',linewidth,'MarkerSize',markersize);
  hold off 
  
 set(handaxes2, 'Box', 'off')
 grid on
 xlim([260 300])
  ylim([60.5 64.5])
  
  handaxes3 = axes('position',[0.42,0.785,0.12,0.12]);
      fontsize = 7;
     linewidth = 0.5;
%  markersize = 0.5;
plot(L, Q_sca,'.-','color',[1.0 0.564 0.141],'LineWidth',linewidth,'MarkerSize',markersize);
 set(handaxes3, 'Box', 'off')
 grid on
 xlim([200 500])
 
 
 handaxes4 = axes('position',[0.47,0.48,0.12,0.12]);
      fontsize = 7;
     linewidth = 0.5;
%  markersize = 0.5;
plot(L, sigma_sca,'.-','color',[0.7 0.8 0.3],'LineWidth',linewidth,'MarkerSize',markersize);
 set(handaxes4, 'Box', 'off')
 grid on
 xlim([200 400])
%   ylim([60 70])
%   ylim([0.6 0.78])
  
%  I = I_o exp(-(b_sca + b_abs)l)   for number of particles
% b_x = <sigm_x> (particle number/m^3)




elseif display_by_Mie == 5
% Optical index reference:  S. Babar, J. H. Weaver; Optical constants of...
% Cu, Ag, and Au revisited; APPLIED OPTICS, Vol. 54, No. 3,...
% http://dx.doi.org/10.1364/AO.54.000477

op_index  = xlsread('Ref_index_Cu_Ag_Au.xlsx');
energy_eV_for_lambda= op_index(:,1);

Wavelength = 1.2398393589807376948907114057097e-6 *10^9 ./ energy_eV_for_lambda; %*10^9 for gettingin nm
Ag_epsi1= op_index(:,6);
Ag_epsi2= op_index(:,7);
Ag_reflectance= op_index(:,5);

interploation_Points = 10^6;

wavelength = zeros(interploation_Points, 1);
Ag_ep1 = zeros(interploation_Points, 1);
Ag_ep2 = zeros(interploation_Points, 1);
Ag_R = zeros(interploation_Points, 1);

    Wavelength(70:end, :) = [];
    wv = Wavelength;
    Agep1 = op_index(:,6);
    Agep2 = op_index(:,7);
    AgR = op_index(:,5);
  wv(isinf(wv)) = nan;
  fillmissing(wv, 'linear'); 
  wv = fillmissing(wv, 'linear');


   X = linspace(min(wv),max(wv),interploation_Points);
   [C,~,idwv] = unique(wv,'stable');
   val1 = accumarray(idwv,Agep1,[],@mean);
   val2 = accumarray(idwv,Agep2,[],@mean);
   val3 = accumarray(idwv,AgR,[],@mean);
   Y = interp1(C,val1,X,'linear','extrap'); 
   Z = interp1(C,val2,X,'linear','extrap');
   V = interp1(C,val3,X,'linear','extrap');
% wavelength = wavelength + X';
% Ag_ep1 = Ag_ep1 + Y';
wavelength = X;
Ag_ep1 = Y;
Ag_ep2 = Z;  
Ag_R = V;

% epsi = epsi1 +i epsi2; m = n+ik => m^2 =(n+ik)^2=n^2+k^2-i2nk = epsil1+ i(epsil2)
% here, abs(epsl)= sqrt((epsil1)^2+(epsil2)^2)
% n = sqrt((abs(epsl)+(epsil1))/2) ; k = sqrt((abs(epsl)-(epsil1))/2)
ep_abs =  sqrt((Ag_ep1).^2+(Ag_ep2).^2);
real_n = sqrt(((ep_abs)+(Ag_ep1))./2);
kappa = sqrt(((ep_abs)-(Ag_ep1))./2);


M_p = (real_n + kappa.*1i);
m_m = 1; % optical index of medium
m_p = 1/m_m .*(linspace(min(M_p),max(M_p),interploation_Points));
%   m_p = (0.135 + 3.99i);

   Lambda = linspace(min(wavelength),max(wavelength),interploation_Points);
%    lambda =200: 1: 1000;
   L= Lambda;% in nm
for KK = 1 : length(L)
  this_L = L(KK);
  this_m_p = m_p(KK);

k = (2* pi*m_m ./this_L);
 a = 5;   % nm
% a = 0.2725846710334890; 
x = (a .* k);
z = (a.*this_m_p .* 2.*pi./this_L);

  n_max = round(x + 4.05 .*(x) .^(1/3) + 2);
 n = 1:n_max; 
%spherical bessel (n,z) = besselj(n+1/2,z)*sqrt(pi/(2*z))
%shi_n(x) = rho * j_n(x)
prefacx = x .*(pi./ (2*x)).^(1/2);
prefacz = z .*(pi./ (2*z)).^(1/2);
shi_x = prefacx .* besselj(n+0.5,x);
kai_x = prefacx .* (besselj(n+0.5,x)+1i*bessely(n+ 0.5,x));
shi_z = prefacz .* besselj(n+0.5,z);
y0 = prefacx .* bessely (n + 0.5, x);
y1 = [-cos(x), y0(1 : n_max-1)];

shi_xn = [sin(x), shi_x(1 : n_max-1)];
shi_zn = [sin(z), shi_z(1 : n_max-1)];

dshi_x = [(shi_xn - n/x.*shi_x)] ;
%  d(kai) =  kai(n-1,x)- n*kai(n,x)/x. ;
dkai_x =  (shi_xn + 1i * y1)- n./x.*(shi_x + 1i*y0);
dshi_z =  (shi_zn - n/z.*shi_z);

 a_n = (this_m_p .* shi_z .* dshi_x - shi_x .* dshi_z)./(this_m_p .* shi_z .* dkai_x - kai_x .* dshi_z);
 b_n = (shi_z .* dshi_x - this_m_p .* shi_x .* dshi_z)./(shi_z .* dkai_x - this_m_p .* kai_x .* dshi_z);

 
sigma_sca(KK) = 2*pi./ k .^2 .* (sum ((2*n+1) .* (abs(a_n) .* abs(a_n)  +  abs(b_n) .* abs(b_n))));
sigma_ext(KK) = 2*pi./ k .^2 .* sum ((2*n+1) .* real(a_n + b_n));
sigma_abso(KK) = (sigma_ext(KK) - sigma_sca(KK));
Q_sca(KK) = sigma_sca(KK)/(pi * a^2);
Q_ext(KK) = sigma_ext(KK)/(pi * a^2);
Q_abso(KK) = sigma_abso(KK)/(pi * a^2);
 
this_L = this_L + 1;
this_m_p=this_m_p+1;
KK = KK+1;

ZP  = [shi_z length(n_max)./shi_z length(n_max -1)];
ZA = dshi_z ./shi_z ;  
end
 
fontsize = 14;
linewidth = 1.5;
 markersize = 1;
 
  subplot(2,1,1)
plot(L, Q_sca,'.-','color',[1.0 0.564 0.141],'LineWidth',linewidth,'MarkerSize',markersize);
  hold on
plot(L, Q_ext,'c-','LineWidth',linewidth,'MarkerSize',markersize);
plot(L, Q_abso,'.-','color',[0.0 0.4 0.0],'LineWidth',linewidth,'MarkerSize',markersize);
  hold off 
  legend('Q_{sca}','Q_{ext}','Q_{abso}','All','Location','northeast');
grid on
    set(gca,'FontSize',fontsize)
    set(gcf,'Color','w')
    ylabel('Q_{sca},Q_{ext},Q_{abso} ')
    xlabel('\lambda (nm)')  
    xlim([200 1000])
    
     
 
 subplot(2,1,2)
plot(L, sigma_sca,'.-','color',[0.7 0.8 0.3],'LineWidth',linewidth,'MarkerSize',markersize);
  hold on
plot(L, sigma_ext,'.-','color',[0.6 0.1 1.0],'LineWidth',linewidth,'MarkerSize',markersize);
plot(L, sigma_abso,'.-','color',[0.0 0.4 0.5],'LineWidth',linewidth,'MarkerSize',markersize);
  hold off 
  legend('\sigma_{sca}','\sigma_{ext}','\sigma_{abso}','All','Location','northeast');
grid on
    set(gca,'FontSize',fontsize)
    set(gcf,'Color','w')
    ylabel('\sigma_{sca},\sigma_{ext},\sigma_{abso} in mm^2')
    xlabel('\lambda (nm)')     
 annotation('textbox',[0.34 0.939 0.7 0.024],'String',{'Mie scattering for Silver NP radius  5 nm of optical index as a function of \lambda '},'FitBoxToText','on');   
  xlim([200 1000])
  
  
  handaxes1 = axes('position',[0.62,0.625,0.2,0.27]);
      fontsize = 7;
     linewidth = 0.5;
%  markersize = 0.5;
plot(L, Q_sca,'.-','color',[1.0 0.564 0.141],'LineWidth',linewidth,'MarkerSize',markersize);
  hold on
plot(L, Q_ext,'c-','LineWidth',linewidth,'MarkerSize',markersize);
plot(L, Q_abso,'.-','color',[0.0 0.4 0.0],'LineWidth',linewidth,'MarkerSize',markersize);
  hold off 
  set(handaxes1, 'Box', 'off')
 grid on
 
xlim([300 330])
 ylim([0.93 0.97])
 
 
  handaxes2 = axes('position',[0.62,0.151,0.2,0.26]);
      fontsize = 7;
     linewidth = 0.5;
%  markersize = 0.5;
plot(L, sigma_sca,'.-','color',[0.7 0.8 0.3],'LineWidth',linewidth,'MarkerSize',markersize);
  hold on
plot(L, sigma_ext,'.-','color',[0.6 0.1 1.0],'LineWidth',linewidth,'MarkerSize',markersize);
plot(L, sigma_abso,'.-','color',[0.0 0.4 0.5],'LineWidth',linewidth,'MarkerSize',markersize);
  hold off 
  set(handaxes2, 'Box', 'off')
 grid on
 xlim([300 330])
 ylim([72 76.5])
 
 
 handaxes3= axes('position',[0.44,0.655,0.15,0.15]);
      fontsize = 7;
     linewidth = 0.5;
%  markersize = 0.5;
plot(L, Q_sca,'.-','color',[1.0 0.564 0.141],'LineWidth',linewidth,'MarkerSize',markersize);
set(handaxes3, 'Box', 'off')
 grid on
 xlim([200 500])
%  ylim([0.93 0.97])

handaxes4= axes('position',[0.45,0.185,0.15,0.15]);
      fontsize = 7;
     linewidth = 0.5;
%  markersize = 0.5;
plot(L, sigma_sca,'.-','color',[0.7 0.8 0.3],'LineWidth',linewidth,'MarkerSize',markersize);
  set(handaxes4, 'Box', 'off')
 grid on
 xlim([200 500])
%  ylim([72 76.5])


%  I = I_o exp(-(b_sca + b_abs)l)   for number of particles
% b_x = <sigm_x> (particle number/m^3)


end