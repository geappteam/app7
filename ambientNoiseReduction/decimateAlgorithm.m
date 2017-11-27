function [ decimatedSample ] = decimateAlgorithm( sample, L )
%decimateAlgorithm : Decimate function of Matlab but recoded
    decimatedSample = zeros(length(sample)/L,1);
    for i = 1:length(sample)/L
        decimatedSample(i) = decimatedSample(i) + sample(L*i-(L-1));
    end
end

