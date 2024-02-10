%buğrahan serttaş
%kanal durum tahmini hesaplaması için yapılır
%SSS değerlerini bulduktan sonra yapılır

function ssshest=chestSSS(enb,rxgrid,startSubframe,sssInd,sssSym0,sssSym5)
dims=lteDLResourceGridSize(enb,1);
L=dims(2);
nSubframes=size(rxgrid,2)/L;
N=ceil(nSubframes/5);
R=size(rxgrid,3);
ssshest=zeros(N,R,numel(sssInd));
% For two phases, one for subframe 0 and one for subframe 5
for phase=1:2
    % Determine which subframe (0 or 5) corresponds to which phase (1
    % or 2) depending on 'startSubframe' (0 or 5).
    if (xor(startSubframe==0,phase==2))
        sym=sssSym0;
    else
        sym=sssSym5;
    end
    % Run all correlations for this phase; correlations are indexed
    % into 'ssshest' with 'idx' and subframes are extracted from the
    % received grid using 'nsf' which is calculated from 'idx'.
    for idx=phase:2:N
        nsf=(idx-1)*5;
        for r=1:R
            subframe=rxgrid(:,(nsf*L)+(1:L),r);
            ssshest(idx,r,:)=subframe(sssInd).*conj(sym);
        end
    end
end

end