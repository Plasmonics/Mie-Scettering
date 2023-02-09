% Mie scattering
%% ------------------------------------
% Angle dependent Mie scattering polar plot
%  Last update: 10/14/2018
%  Author:Hosna Sultana

%% ------------------------------------

theta = 0:1:360; 
%   theta = 40;

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

 for KK = 1 : length(theta)
    
   this_theta = theta(KK);
 p=0; p(1)=1; p(2) = cosd(this_theta);
 t=0; t(1)= cosd(this_theta); t(2)= 3*cosd(2*this_theta);

for i=3:11
    
 pi1 = ((2*i-1)/(i-1))*cosd(this_theta).*p(i-1)-(i/(i-1)).* p(i-2);
% p(i+1) = pi1;
p(i) = pi1;
tau1 = i* cosd(this_theta).*p(i)- (i+1).*p(i-1);
% t(i+1) = tau1;
t(i) = tau1;
end

s11(KK)= sum ((2.* n+1).*(abs(a_n .* p + b_n .* t) .* abs(a_n .* p + b_n .* t))./(n .*(n+1)));
s22(KK)= sum ((2.* n+1).*(abs(a_n .* t + b_n .* p) .* abs(a_n .* t + b_n .* p))./(n .*(n+1)));
 
KK = KK+1;
 end
surS11 = (abs(s11)).* (abs(s11)) ;
surS22 = (abs(s22)).* (abs(s22)) ;
s11_up = (s11(:,(1:181))); theta_up = (theta(:,(1:181)));
s22_down = (s22(:,(181:361))); theta_down = (theta(:,(181:361)));

Sigma_sca_theta= (s11+s22).*(lambda^2/ 8* pi^2);


 ang_up = 0:pi/180:pi ;   ang_down = pi:pi/180:2*pi;  
    
polar(ang_up,s11_up,'c-'); 
 hold on
 polar(ang_down,s22_down,'m-');
 hold off
legend('s11(\theta)', 's22(\theta)','All');
grid on
% xlabel('\theta in degree'); ylabel('s_{11}(\theta) and s_{22}(\theta)')
    
   
% plot( ang,Sigma_sca_theta)
% ZP  = [shi_z length(n_max)./shi_z length(n_max -1)];
% ZA = dshi_z ./shi_z ; 