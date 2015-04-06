close all;
clear all;
clc;

N = 8192;
x = [];
for k = 0:1:N-1
    x(k+1) = k*2 + (k*2+1)*i;
end

y = fft(x);
y = y.';


r = load('test.log');
r = r(:,1) + r(:,2)*i;
max_real = max(abs(real(r)));
max_imag = max(abs(imag(r)));

figure;
plot(abs(real(r) - real(y))./max_real);
grid;

figure;
plot(abs(imag(r) - imag(y))./max_imag);
grid;