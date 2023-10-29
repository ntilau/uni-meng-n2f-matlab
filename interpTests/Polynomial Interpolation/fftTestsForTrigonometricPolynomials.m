%%% fft tests
clc; clear all;
N = 30; % y(x) sampling points
fact = 5; % integer of over-reconstruction M = fact * N
freq = 5;
truncLevel = 1e-3; % truncation level in coefficients preserved
%% tested function
x=linspace(0,2*pi*(N-1)/N,N);
% y=cos(x*freq);
y=exp(1i*x*freq).*cos(x);
figure(1); subplot(2,2,1); 
plot(x,real(y),'.-',x,imag(y),'.-'); axis('tight');
title('Original function');
%% fft function call
Y=fft(y);
X=-floor(N/2):floor(N/2)-1;
subplot(2,2,2); plot(X,real(fftshift(Y)),'.-',X, imag(fftshift(Y)),'.-'); axis('tight');
title('Spectrum');
%% implemented ifft with possible truncation
M=fact*N;
x_rec=linspace(0,2*pi*(M-1)/M,M);
y_rec = zeros(size(x_rec));
angles = x_rec;
% positive freq coeffs
for i=0:floor(N/2)+1
   if (abs(Y(i+1)/max(Y)) > truncLevel) 
      y_rec = y_rec + 1/N*Y(i+1)*exp(1i*angles*i);
   end
end
% negative freq coeffs
for i=floor(N/2)+2:N
   if (abs(Y(i)/max(Y)) > truncLevel)
      y_rec = y_rec + 1/N*Y(i)*exp(-1i*angles*(N+1-i));
   end
end
subplot(2,2,[3 4]); plot(x_rec,real(y_rec),'.-',x_rec,imag(y_rec),'.-');
axis('tight'); title(['Reconstructed function - error : ', ...
   num2str(norm(y-y_rec(1:fact:end))/norm(y)),' with truncation ', ...
   num2str(truncLevel)]);