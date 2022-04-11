%?? IL + I3 + In + d(V3-0)/dt)*Cn = 0  *********

%11 IR3 = -V3*G3 --> IR3 + V3*G3 = 0




% i) Time Domain
%1 V1 = Vin;
%2 (d(V1-V2)/dt)*C1 + (V1-V2)*G1 + Iin = 0
%3 (d(V2-V1)/dt)*C1 + (V2-V1)*G1 + V2*G2 + IL = 0
%4 V2-V3 = L*(d(IL)/dt) --> V2-V3 - 1i*omega*(L)*IL = 0
%5 I3 = V3*G3
%6 IL+I3 = 0
%7 V4 = a*I3
%8 (V4-Vo)*G4 + I4 = 0
%9 (Vo-V4)*G4 + Vo*Go = 0
%10 I3 = In + V3*G3 + d(V3-0)/dt)*Cn --> In + d(V3-0)/dt)*Cn + V3*G3 - I3 = 0


% G(10,3) = G3; G(10,8) = -1; G(10,10) = 1;

% Pin = abs(FFTin/n).^2;
% Pout = abs(FFTout/n).^2;




% i) Time Domain
% 1 V1 = Vin;
% 2 (d(V1-V2)/dt)*C + (V1-V2)*G1 + Iin = 0
% 3 (d(V2-V1)/dt)*C + (V2-V1)*G1 + V2*G2 + IL = 0
% 4 V2-V3 = L*(d(IL)/dt)  % V2-V3 - 1i*omega*(L)*IL = 0
% 5 I3 = V3*G3
% 6 IL+I3 = 0
% 7 V4 = a*I3
% 8 (V4-Vo)*G4 + I4 = 0
% 9 (Vo-V4)*G4 + Vo*Go = 0
% 
% 
% ii) Frequency Domain
% 1 V1 = Vin;
% 2 (V1-V2)*1i*omega*C + (V1-V2)*G1 + Iin = 0
% 3 (V2-V1)*1i*omega*C + (V2-V1)*G1 + V2*G2 + IL = 0
% 4 V2-V3 - 1i*omega*(L)*IL = 0
% 5 V3*G3 - I3 = 0
% 6 IL+I3 = 0
% 7 V4 - a*I3 = 0
% 8 (V4-Vo)*G4 + I4 = 0
% 9 (Vo-V4)*G4 + Vo*Go = 0


% i) Time Domain
% 1 V1 = Vin;
% 2 (d(V1-V2)/dt)*C1 + (V1-V2)*G1 + Iin = 0
% 3 (d(V2-V1)/dt)*C1 + (V2-V1)*G1 + V2*G2 + IL = 0
% 4 V2-V3 = L*(d(IL)/dt) --> V2-V3 - 1i*omega*(L)*IL = 0
% 5 I3 = V3*G3
% 6 IL + I3 + In + d(V3-0)/dt)*Cn= 0
% 7 V4 = a*I3
% 8 (V4-Vo)*G4 + I4 = 0
% 9 (Vo-V4)*G4 + Vo*Go = 0


