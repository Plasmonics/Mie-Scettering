% Mie scattering
%% ------------------------------------
% Angle dependent Mie scattering cross-section 
%  Last update: 10/14/2018
%  Author:Hosna Sultana

%% ------------------------------------

ang = 0:1:180;   
 % ang = 0;
 %ang = 0:2*pi/360:2*pi;
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

for KK = 1 : length(ang)
    
  this_ang = ang(KK);

p(1,:) = ones(1,size(this_ang,2));
t(1,:) = cosd(this_ang);
p(2,:) = 3*cosd(this_ang);
t(2,:) = 2*cosd(this_ang).*p(2,:)-3;
for n=3:n_max
p(n,:) = ((2*n-1)*cosd(this_ang).*p(n-1,:) - n*p(n-2,:))/(n-1);
t(n,:) = n*cosd(this_ang).*p(n,:) - (n+1)*p(n-1,:);
end
 %s1(KK)= sum ((2*n+1)*(abs(a_n.*(p)' + b_n.* (t)').* abs(a_n.*(p)' + b_n.* (t)'))./(n*(n+1)));
 %s2(KK)= sum ((2*n+1)*(abs(a_n.*(t)' + b_n.* (p)').* abs(a_n.*(t)' + b_n.* (p)'))./(n*(n+1)));
 %s1(KK)= 1i*( sum ((2*n+1)*(a_n.*(p)'+ b_n.* (t)'))./k.*(n*(n+1)));
%s2(KK)= 1i*(sum ((2*n+1)*(a_n.*(t)' + b_n.* (p)')./k.*(n*(n+1))));

 s1(KK)= ( sum ((2*n+1)*(a_n.*(p)'+ b_n.* (t)'))./(n*(n+1)));
s2(KK)= (sum ((2*n+1)*(a_n.*(t)' + b_n.* (p)')./(n*(n+1))));

 sq_s1 = abs((s1)).^2 ;
sq_s2 = abs((s2)).^2 ;
 sq_s11 =0.5 .* (abs((s1)).^2 + abs(s2).^2);
sq_s12 = 0.5 .* (abs((s1)).^2 - abs(s2).^2); 

% S1= sum (((2*n+1)/(n*(n+1))*(abs(a_n.*p + b_n.* t).* abs(a_n.*p + b_n.* t))));

% S2= sum (((2*n+1)/(n*(n+1)))*(abs(a_n.*t + b_n.* p).* abs(a_n.*t + b_n.* p)));

%  Sigma_sca_theta(KK)= (s11+s22).*(lambda^2/ 8* pi^2);

 this_ang = this_ang + 1;
KK = KK+1;
end
surS11 = (abs(s1)).* (abs(s1)) ;
surS22 = (abs(s2)).* (abs(s2)) ;
s11= (1/2)* (surS11 + surS22);
s12= (1/2)* (surS22 - surS11);


                           display_by_plot = 1; 
%                           display_by_plot = 2; 
%                           display_by_plot = 3; 

 
                         if  display_by_plot == 1
  
                         subplot (2,1,1)                   
 plot( ang,sq_s1)
 hold on
 plot( ang,sq_s2)
 hold off
legend('\mid s_{1}^{2}(\theta)\mid',' \mid s_{2}^{2}(\theta)\mid','All','Location','northeast');
ylabel('\mid s_{1}^{2}(\theta)\mid and \mid s_{2}^{2}(\theta)\mid');
subplot(2,1,2)
plot( ang,sq_s11)
 hold on
 plot( ang,sq_s12)
 hold off
legend('S_{11}(\theta)','S_{12}(\theta)','All','Location','northeast');
xlabel('\theta in degree'); ylabel('S_{11}(\theta)and s_{12}(\theta)');

                                 elseif  display_by_plot == 2

 plot( ang,s1)
 hold on
 plot( ang,s2)
 hold off
legend('s_{11}','s_{22}','All','Location','northeast');


grid on
                               elseif display_by_plot ==3


 s11_up = (s11(:,(1:181))); theta_up = (theta(:,(1:181)));
s22_down = (s22(:,(181:361))); theta_down = (theta(:,(181:361)));
 ang_up = 0:pi/180:pi ;   ang_down = pi:pi/180:2*pi;  
    
polar(ang_up,s11_up,'c-'); 
 hold on
 polar(ang_down,s22_down,'m-');
 hold off
legend('s11(\theta)', 's22(\theta)','All');
grid on
% xlabel('\theta in degree'); ylabel('s_{11}(\theta) and s_{22}(\theta)')
    
 end     