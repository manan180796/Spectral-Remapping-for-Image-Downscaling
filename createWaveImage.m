function [Result] = createWaveImage(wave,imageSize)
if(wave(2)<0)
    wave(2)=wave(3)+1;
end
if(wave(3)<0)
    wave(3)=wave(3)+1;
end
a = int16(imageSize(2)*wave(2));
b = int16(imageSize(1)*wave(3));
waveFFT = zeros(imageSize);
waveFFT(b+1,a+1) = wave(1)*imageSize(1)*imageSize(2);
waveFFT(b+1,a+1)=waveFFT(b+1,a+1)*exp(2i*pi*wave(4));
Result = ifft2(waveFFT);
end