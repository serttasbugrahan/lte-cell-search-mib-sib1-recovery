%buğrahan serttaş
%bu işlemler ana kodun ilk başındaki RB değerine göre oluşan kodlardır
ofdmInfo.SamplingRate = '';
ofdmInfo.Windowing = '';
ofdmInfo.CyclicPrefixLengths = '';
ofdmInfo.Nfft = '';

if enb.NDLRB == 6
    ofdmInfo.SamplingRate = 1920000;
elseif enb.NDLRB == 15
    ofdmInfo.SamplingRate = 3840000;
elseif enb.NDLRB == 25
    ofdmInfo.SamplingRate = 7680000;
elseif enb.NDLRB == 50
    ofdmInfo.SamplingRate = 15360000;
elseif enb.NDLRB == 75
    ofdmInfo.SamplingRate = 30720000;
elseif enb.NDLRB == 100
    ofdmInfo.SamplingRate = 30720000;
end

if enb.NDLRB == 6
    ofdmInfo.Nfft = 128;
    ofdmInfo.Windowing = 4;
elseif enb.NDLRB == 15
    ofdmInfo.Nfft = 256;
    ofdmInfo.Windowing = 6;
elseif enb.NDLRB == 25
    ofdmInfo.Nfft = 512;
    ofdmInfo.Windowing = 4;
elseif enb.NDLRB == 50
    ofdmInfo.Nfft = 1024;
    ofdmInfo.Windowing = 6;
elseif enb.NDLRB == 75
    ofdmInfo.Nfft = 2048;
    ofdmInfo.Windowing = 8;
elseif enb.NDLRB == 100
    ofdmInfo.Nfft = 2048;
    ofdmInfo.Windowing = 8;
end


if ofdmInfo.Nfft == 128  && strcmp(enb.CyclicPrefix, 'Normal')
    ofdmInfo.CyclicPrefixLengths = [10 9 9 9 9 9 9 10 9 9 9 9 9 9];
elseif ofdmInfo.Nfft == 128  && strcmp(enb.CyclicPrefix, 'Extended')
    ofdmInfo.CyclicPrefixLengths = [32 32 32 32 32 32 32 32 32 32 32 32];
elseif ofdmInfo.Nfft == 256  && strcmp(enb.CyclicPrefix, 'Normal')
    ofdmInfo.CyclicPrefixLengths = [20 18 18 18 18 18 18 20 18 18 18 18 18 18];
elseif ofdmInfo.Nfft == 256  && strcmp(enb.CyclicPrefix, 'Extended')
    ofdmInfo.CyclicPrefixLengths = [64 64 64 64 64 64 64 64 64 64 64 64];
elseif ofdmInfo.Nfft == 512  && strcmp(enb.CyclicPrefix, 'Normal')
    ofdmInfo.CyclicPrefixLengths = 	[40 36 36 36 36 36 36 40 36 36 36 36 36 36];
elseif ofdmInfo.Nfft == 512  && strcmp(enb.CyclicPrefix, 'Extended')
    ofdmInfo.CyclicPrefixLengths = [128 128 128 128 128 128 128 128 128 128 128 128];
elseif ofdmInfo.Nfft == 1024  && strcmp(enb.CyclicPrefix, 'Normal')
    ofdmInfo.CyclicPrefixLengths = 	[80 72 72 72 72 72 72 80 72 72 72 72 72 72];
elseif ofdmInfo.Nfft == 1024 && strcmp(enb.CyclicPrefix, 'Extended')
    ofdmInfo.CyclicPrefixLengths = [256 256 256 256 256 256 256 256 256 256 256 256];
elseif ofdmInfo.Nfft == 2048  && strcmp(enb.CyclicPrefix, 'Normal')
    ofdmInfo.CyclicPrefixLengths = 	[160 144 144 144 144 144 144 160 144 144 144 144 144 144];
elseif ofdmInfo.Nfft == 2048 && strcmp(enb.CyclicPrefix, 'Extended')
    ofdmInfo.CyclicPrefixLengths = [512 512 512 512 512 512 512 512 512 512 512 512];
end