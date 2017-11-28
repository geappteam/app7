function [ output_args ] = customFreqz( b, a, minW, maxW)
    
    nbPoints = 1000;

    if ~exist('maxW','var')
        maxW = 1;
    end
    if ~exist('minW','var')
        minW = 0;
    end
    
    incW = (maxW - minW) / (nbPoints);
    
    W = (minW:incW:maxW-incW);
    teta = W*pi;
    
    z = exp(-1j*teta);
    
    rFreq = polyval(b,z)./polyval(a,z);
    
    figure()
    subplot(2,1,1)
    plot(W, mag2db(abs(rFreq)))
    title('Reponse en frequence (freqz maison)')
    xlabel('Normalized Frequency (* \pi rad/sample)')
    ylabel('Magnitude (dB)')
    subplot(2,1,2)
    plot(W, -rad2deg(angle(rFreq)))
    xlabel('Normalized Frequency (* \pi rad/sample)')
    ylabel('Phase (degrees)')
    output_args = rFreq;
end

