%buğrahan serttaş
%subframe 0'daki SSS değerleri içindir
%SSS indis değerlerinden sonra gerçekleşir

function SSS_subfram0 = LTE_SSS_0(enb)
slot_num = 0;
x = enb.NCellID;
n_id_1 = [];
n_id_2 = [];

for b = 1:2 % b değerleri aralığı için bir döngü başlatın
    a = (x - b) / 3;
    if mod(a, 1) == 0 &&  mod(b, 1) == 0 % eğer a ve b tam sayıysa
        n_id_1 = [n_id_1 a]; % a değerlerini sakla
        n_id_2 = [n_id_2 b]; % b değerlerini sakla
    end
end
qp=floor(n_id_1/30);
q=floor((n_id_1+qp*(qp+1)/2)/30);
mp=n_id_1+q*(q+1)/2;
m0=mod(mp,31);
m1=mod(m0+floor(mp/31)+1,31);

%s_td=[0 0 0 0 1];
%for t=1:26
%  s_td=[s_td mod(s_td(end-2)+s_td(end-4),2)];
%end
s_td=[0 0 0 0 1 0 0 1 0 1 1 0 0 1 1 1 1 1 0 0 0 1 1 0 1 1 1 0 1 0 1];
s_td=1-2*s_td;

%c_td=[0 0 0 0 1];
%for t=1:26
%  c_td=[c_td mod(c_td(end-1)+c_td(end-4),2)];
%end
c_td=[0 0 0 0 1 0 1 0 1 1 1 0 1 1 0 0 0 1 1 1 1 1 0 0 1 1 0 1 0 0 1];
c_td=1-2*c_td;

%z_td=[0 0 0 0 1];
%for t=1:26
%  z_td=[z_td mod(z_td(end)+z_td(end-2)+z_td(end-3)+z_td(end-4),2)];
%end
z_td=[0 0 0 0 1 1 1 0 0 1 1 0 1 1 1 1 1 0 1 0 0 0 1 0 0 1 0 1 0 1 1];
z_td=1-2*z_td;

s0_m0=s_td(mod(m0:30+m0,31)+1);
s1_m1=s_td(mod(m1:30+m1,31)+1);

c0=c_td(mod(n_id_2:30+n_id_2,31)+1);
c1=c_td(mod(n_id_2+3:30+n_id_2+3,31)+1);

z1_m0=z_td(mod((0:30)+mod(m0,8),31)+1);
z1_m1=z_td(mod((0:30)+mod(m1,8),31)+1);

if (slot_num==0)
    SSS_subfram0(2:2:62)=s1_m1.*c1.*z1_m0;
    SSS_subfram0(1:2:62)=s0_m0.*c0;
    SSS_subfram0 = transpose(SSS_subfram0);
    
elseif (slot_num==10)
    SSS_subfram5(2:2:62)=s0_m0.*c1.*z1_m1;
    SSS_subfram5(1:2:62)=s1_m1.*c0;
    SSS_subfram5 = transpose(SSS_subfram5);
end
end