%buğrahan serttaş
%PCFICH indis değerlerini bulmak için yapılır
%PCFICH işleminin başında kullanılır

function pcfich_Indices = lte_PCFICH_Indices(enb)
N_RB_SC=12;
pcfich_Indices = uint32(zeros(16, 2));
k = (N_RB_SC/2) * mod(enb.NCellID, 2*enb.NDLRB);
k1= k + round(enb.NDLRB/2)*N_RB_SC/2;
k2=k + 2*round(enb.NDLRB/2)*N_RB_SC/2;
k3=k + 3*round(enb.NDLRB/2)*N_RB_SC/2;
pcfich_Indices(1:16,1) = [k+1,k+2,k+4,k+5,k1+1,k1+2,k1+4,k1+5,k2+1,k2+2,k2+4,k2+5,k3+1,k3+2,k3+4,k3+5];
pcfich_Indices(1:16,2) = [k+1+8400,k+2+8400,k+4+8400,k+5+8400,k1+1+8400,k1+2+8400,k1+4+8400,k1+5+8400,k2+1+8400,k2+2+8400,k2+4+8400,k2+5+8400,k3+1+8400,k3+2+8400,k3+4+8400,k3+5+8400];

end