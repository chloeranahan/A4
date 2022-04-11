clear

R1 = 1;
C1 = 0.25;
R2 = 2;
L = 0.2;
R3 = 33.29; %taken from linear fit from Q1 & Q2
a = 100;
R4 = 0.1;
Ro = 1000;

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
Go = 1/Ro;


% i) Time Domain
%1 V1 = Vin;
%2 (d(V1-V2)/dt)*C + (V1-V2)*G1 + Iin = 0
%3 (d(V2-V1)/dt)*C + (V2-V1)*G1 + V2*G2 + IL = 0
%4 V2-V3 = L*(d(IL)/dt)  % V2-V3 - 1i*omega*(L)*IL = 0
%5 I3 = V3*G3
%6 IL+I3 = 0
%7 V4 = a*I3
%8 (V4-Vo)*G4 + I4 = 0
%9 (Vo-V4)*G4 + Vo*Go = 0



%iii) Matrices
omega = 0;
N = 9;

%X = [V1 V2 V3 V4 Vo Iin IL I3 I4];
n3 = 3;
no = 5;

G = zeros(N);  
F = zeros(N,1);
C = zeros(N);
CFinal = C.*(1i*omega);


G(1,1) = 1;

G(2,1) = G1; G(2,2) = -G1; G(2,6) = 1;
C(2,1) = C1; C(2,2) = -C1;

G(3,1) = -G1; G(3,2) = G1 + G2; G(3,7) = 1;
C(3,1) = -C1; C(3,2) = C1;

G(4,2) = 1; G(4,3) = -1;
% C(4,7) = -L;
C(4,7) = L;

G(5,3) = G3; G(5,8) = -1;

G(6,7) = 1; G(6,8) = 1;

G(7,4) = 1; G(7,8) = -a;

G(8,4) = G4; G(8,5) = -G4; G(8,9) = 1;

G(9,4) = -G4; G(9,5) = G4 + Go;


Vin = [];
Vo = [];
Vp = zeros(N,1);

tstop = 1;
tstep = 1/1000;

inputA = 0;
inputB = 0;
inputC = 1;

freq = 1/0.03;
% freq = 1/0.3;
% freq = 1/0.009;
% freq = 1/0.09;

std = 0.03;
delay = 0.06;

%x = linspace(0,1000,1000);
mu = 560;
sig = 30;
amp = 1;
vo = 0;
gaus = @(i,mu,sig,amp,vo)amp*exp(-(((i-mu).^2)/(2*sig.^2)))+vo;

T = linspace(0,tstop,1000);

for i = 1:1000
    t(i) = i*tstep;
    
    if t(i) < 0.03 && inputA == 1
        Vin(i) = 0;
    elseif t(i) >= 0.03 && inputA == 1
        Vin(i) = 1;
        
    elseif inputB ==1
        Vin(i) = sin(2*pi*freq*t(i));
        
    elseif inputC == 1
        Vin(i) = gaus(i,mu,sig,amp,vo);
    end
    
    F(1) = Vin(i);
    
    H = C/tstep + G;
    V = H\(F + ((C/tstep)*Vp));
    
    Vo(i) = V(no);
    Vp = V;
end


% iii)
figure(1)
plot(T,Vin)
hold on
plot(T,Vo)
hold off
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('Vin','Vo')
% title('Frequency = 1/(0.09) Hz')


% iv)
% source of sample code: https://www.mathworks.com/help/matlab/ref/fft.html

Fs = 1000;
t_v= T; 
L = length(t_v); % Signal length

n = L;

FFTin = fft(Vin,n);
FFTout = fft(Vo,n);

f =Fs*([(((-n/2)+1):0)/n  (1:(n/2))/n]-1/Fs);

figure(2)
plot(f,fftshift(abs(FFTin)));
hold on
plot(f,fftshift(abs(FFTout))) 
hold off
title('Input and Output Signals in Frequency Domain')
xlabel('Frequency (Hz)')
ylabel('Fourrier Transform (arb. units)')
legend('Input Signal','Output Signal')





% v)
tstep = 1/500;
% tstep = 1/100;
% tstep = 1/50;

Vin = [];
Vo = [];
Vp = zeros(N,1);

tstop = 1;

freq = 1/0.03;

std = 0.03;
delay = 0.06;


mu = 0.56/tstep;
sig = 0.03/tstep;
amp = 1;
vo = 0;
gaus = @(i,mu,sig,amp,vo)amp*exp(-(((i-mu).^2)/(2*sig.^2)))+vo;

T = linspace(0,tstop,1/tstep);

if tstep == 1/500
    for i = 1:500
        t(i) = i*tstep;

        if t(i) < 0.03 && inputA == 1
            Vin(i) = 0;
        elseif t(i) >= 0.03 && inputA == 1
            Vin(i) = 1;

        elseif inputB ==1
            Vin(i) = sin(2*pi*freq*t(i));

        elseif inputC == 1
            Vin(i) = gaus(i,mu,sig,amp,vo);
        end

        F(1) = Vin(i);

        H = C/tstep + G;
        V = H\(F + ((C/tstep)*Vp));

        Vo(i) = V(no);
        Vp = V;
    end
    
elseif tstep == 1/100
    for i = 1:100
    t(i) = i*tstep;
    
    if t(i) < 0.03 && inputA == 1
        Vin(i) = 0;
    elseif t(i) >= 0.03 && inputA == 1
        Vin(i) = 1;
        
    elseif inputB ==1
        Vin(i) = sin(2*pi*freq*t(i));
        
    elseif inputC == 1
        Vin(i) = gaus(i,mu,sig,amp,vo);
    end
    
    F(1) = Vin(i);
    
    H = C/tstep + G;
    V = H\(F + ((C/tstep)*Vp));
    
    Vo(i) = V(no);
    Vp = V;
    end
    
elseif tstep == 1/50
    for i = 1:50
    t(i) = i*tstep;
    
    if t(i) < 0.03 && inputA == 1
        Vin(i) = 0;
    elseif t(i) >= 0.03 && inputA == 1
        Vin(i) = 1;
        
    elseif inputB ==1
        Vin(i) = sin(2*pi*freq*t(i));
        
    elseif inputC == 1
        Vin(i) = gaus(i,mu,sig,amp,vo);
    end
    
    F(1) = Vin(i);
    
    H = C/tstep + G;
    V = H\(F + ((C/tstep)*Vp));
    
    Vo(i) = V(no);
    Vp = V;
end

end


figure(3)
plot(T,Vin)
hold on
plot(T,Vo)
hold off
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('Vin','Vo')
title('500 timesteps')



Fs = 1000;
t_v= T; 
L = length(t_v); % Signal length

n = L;

FFTin = fft(Vin,n);
FFTout = fft(Vo,n);

f =Fs*([(((-n/2)+1):0)/n  (1:(n/2))/n]-1/Fs);

figure(4)
plot(f,fftshift(abs(FFTin)));
hold on
plot(f,fftshift(abs(FFTout))) 
hold off
title('Input and Output Signals in Frequency Domain')
xlabel('Frequency (Hz)')
ylabel('Fourrier Transform (arb. units)')
legend('Input Signal','Output Signal')

