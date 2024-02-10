%buğrahan serttaş
%PSS değerlerini bulmak için yapılır PSS indis işleminden sonra gerçekleşir
function zadoff_chu = lte_PSS_fn(enb)
u_shift = [25 29 34];
% Generate PSS for NID = 0
zadoff_chu = [];
for n = 0:61
    u = u_shift(enb.NCellID+1);
    if n <= 30
        d = exp(-j*pi*u*n*(n+1)/63);
    else
        d = exp(-j*pi*u*(n+1)*(n+2)/63);
    end
    zadoff_chu = [zadoff_chu d];
end
zadoff_chu = zadoff_chu';
end