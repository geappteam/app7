fc = 2000;
Fs = 16000;
wc = fc*2*pi;

% Filter to transfer to z domain
cDen = [1/wc^2, sqrt(2)/wc, 1];
fonctionEnS = tf(1,  cDen)
figure()
bode(fonctionEnS)
title('Analog filter')

% Expected result after transfer as defined by requirements
Wn = fc / (Fs /2);
[expect_b, expect_a] = butter(2, Wn)
figure()
soziFreqz(expect_b,expect_a);
title('Expected filter from butter()');

%% Biliear transform
T = 1/Fs;
cT = 2/T;
C = cDen(1) * cT^2 + cDen(2) * cT + 1;

b = [1, 2, 1] / C
a = [ 1 ,  (2 - cDen(1)* 2 * cT^2) / C,  (cDen(1) * cT^2 - cDen(2) * cT + 1)/ C]

figure()
soziFreqz(b,a);
title('Filter obtained with bilinear transform');