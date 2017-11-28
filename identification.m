% In this script, we identify the effect of a room on a sound signal
% and we apply it onto someone`s voice

% Maximum Impulse response length (To be adjusted)
maxM = 5000;

% Maximum error allowed
maxE = 0.5;

% Sound injected into the room
[pureWN, Fs] = audioread('res/bruit_emis_32kHz.wav');
[pureVoice, ~] = audioread('res/parole_emise_32khz.wav');

% Sound with the effect of the room applied
[taintedWN, ~] = audioread('res/bruit_capte_32kHz.wav');
[taintedVoice, ~] = audioread('res/parole_captee_32khz.wav');



Rxx = xcorr(pureWN, pureWN, maxM);
Rxy = xcorr(pureWN, taintedWN, maxM);

bA = toeplitz(Rxx(maxM+1:maxM*2+1));

E = 0;
m = maxM;
mStep = floor(maxM/2);

while (mStep >= 1)
    
    if (E <= maxE)
        m = m - mStep;
    else
        m = m + mStep;
    end
    
    A = bA(1:m+1,1:m+1);
    b = Rxy(maxM+1:-1:1+maxM-m);
    h = linsolve(A, b);
    
    prediction = filter(h, 1, pureWN);
    E = evalErrorPower(prediction, taintedWN);
    
    if (mStep == 1)
        if (E > maxE)
            continue;
        end
        
        break;
    end
    
    mStep = round(mStep/2);
end

% Outputting the parameters computed
m
N = m+1
E

% Playing the filtered voice sample
out = filter(h, 1, pureVoice);
plyer = audioplayer(out,Fs);
play(plyer);