%buğrahan serttaş
mib = transpose(mib);
mib = char('0' + mib);
if isequal(mib(1,1:3), '000')
    enb.NDLRB = 6;
elseif isequal(mib(1,1:3), '001')
    enb.NDLRB = 15;
elseif isequal(mib(1,1:3), '010')
    enb.NDLRB = 25;
elseif isequal(mib(1, 1:3), '011')
    enb.NDLRB = 50;
elseif isequal(mib(1,1:3), '100')
    enb.NDLRB = 75;
elseif isequal(mib(1,1:3), '101')
    enb.NDLRB = 100;
end
if mib(1,4) == '0'
    enb.PHICHDuration = 'Normal';
elseif mib(1,4) == '1'
    enb.PHICHDuration = 'Extended';
end
if  isequal(mib(1,5:6), '00')
    enb.Ng = 'Sixth';
elseif isequal(mib(1,5:6), '01')
    enb.Ng = 'Half';
elseif isequal(mib(1,5:6), '10')
    enb.Ng = 'One';
elseif isequal(mib(1,5:6), '11')
    enb.Ng = 'Two';
end
enb.NFrame = bin2dec(num2str(mib(7:16)));