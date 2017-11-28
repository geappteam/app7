function [ output_args ] = customFreqz( b, a, minW, maxW)
    
    nbPoints = 1000;

    if ~exist('maxW','var')
        maxW = 1;
    end
    if ~exist('minW','var')
        minW = 0;
    end
    
    incW = (maxW - minW) / (nbPoints);
    
    teta = (minW:incW:maxW)*pi;
    
    z = exp(-1j*teta);
    
    rFreq = polyval(b,z)./polyval(a,z);
    
    figure()
    subplot(2,1,1)
    plot(teta, mag2db(abs(rFreq)))
    title('Reponse en frequence (freqz maison)')
    subplot(2,1,2)
    plot(teta, angle(rFreq))
    
    output_args = rFreq;
end

