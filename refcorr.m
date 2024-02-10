%buğrahan serttaş
%Çapraz korelasyon yapılır(iki sinyal arasındaki benzerlik için önemli birişlemdir)
%PSS işlemleri bittikten sonra yapılır
function corr = refcorr(downsampled,ref)
firstnonzero=find(ref~=0,1);
firstzero=find(ref(firstnonzero:end)==0,1);
corr=flipud(fftfilt(conj(ref(firstnonzero+(0:firstzero-2))),flipud(downsampled(firstnonzero-firstzero+1:end))));
corr=[corr(firstzero:end); zeros(firstnonzero-1,1)];
end