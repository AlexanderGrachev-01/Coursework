clear all;
clc;

tic

a = 10;       
b = 1*a;    
c = 3e+8;
z0 = 2*a;     
d_r = 0.01;
r_int = 0:d_r:b;
phi_int = 0:0.01:2*pi;
Th = (8)*pi/180;
N_FFT = 8192/32; 
T = 25/c;
d_t = T/(N_FFT - 1);с

kx = sin(Th);
ky = 0;
kz = cos(Th);

xs = r_int'*cos(phi_int);
ys = r_int'*sin(phi_int);
zs = z0;
z = zs;   
ro = r_int;  

t_in = z0/c*(1 + kz)/kz - T/2:d_t:z0/c*(1 + kz)/kz + T/2;
t_delay = (kx*xs + kz*zs)/c;
for k1 = 1:length(phi_int)
    E_e1 = zeros(length(r_int), N_FFT);
    t_d = zeros(1, length(r_int));
    for k2 = 1:length(r_int)
        E_e = zeros(1, N_FFT);
        t_d(k2) = round(t_delay(k2, k1)/d_t)*d_t;
         
        t = t_in - t_d(k2);
        B = sqrt((c*t).^2 - z^2);
        if abs(ro(k2)) <= a
           i2 = find((c*t >= z) & (c*t < sqrt(z^2 + (a - abs(ro(k2)))^2)));
           E_e(i2) = ones(1, length(i2));
        end
        i3 = find((c*t >= sqrt(z^2 + (a - abs(ro(k2)))^2)) & (c*t < sqrt(z^2 + (a + abs(ro(k2)))^2)));
        E_e(i3) = 1/pi.*acos((-a^2 + abs(ro(k2))^2 + B(i3).^2)./(2*abs(ro(k2))*B(i3)));
        
        j1 = find(E_e > 0);
        if length(j1) > 0
            E_e1(k2, j1(1)) = E_e(j1(1));
        end
%             plot(t_in*1e+9, E_e1(k2, :));hold on;grid
%             участок(t_in*1e+9, E_e1(k2, :));держит;сетка
    end
%     grid
%     сетка 
    E_b1(k1, :) = trapz(r_int, E_e1.*(r_int'*ones(1, N_FFT)));

%     figure(5)
%     plot(t_in*1e+9, E_b1(k1, :)); hold on; grid
end
% grid
E_b = 1/2/pi/c*trapz(phi_int, E_b1);
E_b_diff = diff(E_b)/d_t;
    
figure(1);
plot((t_in - t_in(1))*1e+9, E_b*1e+9); grid
xlabel('time, nsec');
ylabel('Amplitude');
title('E field recovering. Antiderivative. 2*a = 20 m. \Theta = 3^o');

figure(2);
plot((t_in(1:end - 1) - t_in(1))*1e+9, E_b_diff); grid
xlabel('time, nsec');
ylabel('Amplitude');
title('E field recovering, T_w =  nsec. 2*a = 20 m. z_0 = *a. \Theta = 3^o');

toc



