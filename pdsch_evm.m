%buğrahan serttaş
%pdsch işleminden sonra evm hesaplaması için yapılır
ref_QPSK = refSymbols{1};
%scatterplot(ref_QPSK)
symbols = pdschSymbols{1};
%scatterplot(symbols)

for ii = 1:length(symbols)
    dist(ii) = norm(symbols(ii)-ref_QPSK(ii));
    Adist(ii) = (symbols(ii)-ref_QPSK(ii));
end

peak_EVM = max(dist)*100
rms_EVM = rms(Adist)*100