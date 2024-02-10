%Buğrahan serttaş
%bu işlem PCFICH değerini bulduktan sonra yapılır
if enb.NDLRB == 6 || enb.NDLRB == 15 || enb.NDLRB == 25%bu işlemi ts 36.508 4.3.3.3'e göre yaptım
    cfi = 3;
elseif enb.NDLRB == 50 || enb.NDLRB == 75 || enb.NDLRB == 100
    cfi = 2;
end