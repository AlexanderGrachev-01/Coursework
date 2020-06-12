clear all;
clc;

% Programm of time-space electic fields or radiation pattern calculation
% for circular aperture measured with infinitely small zond

tic
% Design parameters
a = 10;         % Aperture radius, m
b = 3*a;          % Radius of zond surface

lam = a/10;     % Wavelength in free space
c = 3e+8;       % speed of light, m/sec

z0 = 2*a;       % Distance to zond surface

r_int = 0:0.01:b;
phi_int = 0:0.01:2*pi;

Th = (0:0.1:10)*pi/180; % Theta angle
N_Th = length(Th);

N_FFT = 8192/4;  % Number of FFT points (should be varied for faster calculation without accuracy degradation)
T = 200/c*(N_FFT - 1)/N_FFT;    % Time interval for FFT
d_t = T/(N_FFT - 1); % Sample time, sec

xs = r_int'*cos(phi_int);
ys = r_int'*sin(phi_int);
zs = z0;
z = zs;     % Distance from observation point to aperture plane, m  
ro = r_int;   % Vector R projection to aperture plane, m

for i = 1:N_Th
    kx = sin(Th(i));
    ky = 0;
    kz = cos(Th(i));
    t_in = z0/c*(1 + kz)/kz - T/2:d_t:z0/c*(1 + kz)/kz + T/2;
    t_delay = (kx*xs + kz*zs)/c;
    for k1 = 1:length(phi_int)
        E_e1 = zeros(length(r_int), N_FFT);  
        for k2 = 1:length(r_int)
            E_e = zeros(1, N_FFT);
            t_d(k2) = round(t_delay(k2, k1)/d_t)*d_t;
            t = t_in - t_d(k2);
            B = sqrt((c*t).^2 - z^2);
            if abs(ro(k2)) <= a
                i2 = find((c*t >= z) & (c*t < sqrt(z^2 + (a - abs(ro(k2)))^2)));
                E_e(i2) = z^2./((c*t(i2)).^2);
            end
            i3 = find((c*t >= sqrt(z^2 + (a - abs(ro(k2)))^2)) & (c*t < sqrt(z^2 + (a + abs(ro(k2)))^2)));
            E_e(i3) = z^2/pi./((c*t(i3)).^2).*acos((-a^2 + abs(ro(k2))^2 + B(i3).^2)./(2*abs(ro(k2))*B(i3)));
            
            j1 = find(E_e > 0);
            if length(j1) > 0
                E_e1(k2, j1(1)) = E_e(j1(1));
            end
        end
        E_b1(k1, :) = trapz(r_int, E_e1.*(r_int'*ones(1, N_FFT))*cos(Th(i)));
    end
    E_b = 1/2/pi/c*trapz(phi_int, E_b1);
%     
%     E_b_diff = diff(E_b)/d_t;
%     j1 = find(E_b_diff < 0);
%     E_b_diff(j1) = zeros(1, length(j1));
    E_f(i, :) = fft(E_b, N_FFT);
    i
end

if (length(Th) == 1)
    figure(1);
    plot((t)*1e+9, E_b); grid
    xlabel('time, nsec');
    ylabel('Amplitude');
    title('E field recovering. Antiderivative. 2*a = 20 m. \Theta = 3^o');

    figure(2);
    plot((t(1:end - 1))*1e+9, E_b_diff); grid
    xlabel('time, nsec');
    ylabel('Amplitude');
    title('E field recovering. 2*a = 20 m. \Theta = 3^o');
else
    d_f = 1/N_FFT/d_t;
    k = 1:1:N_FFT/2;
    figure(3);
    [F, p] = min(abs(d_f*k - 3e+8/lam));
    plot(2*a/lam*sin(Th), 20*log10(abs(E_f(:, p))/max(abs(E_f(:, p))))); grid
    xlabel('2*a/\lambda*sin(\Theta)');
    ylabel('Amplitude, dB');
    title('Recovered RP of CA for a/\lambda = 10. z0 = 4*a. r = 1.0*a');
end

toc



