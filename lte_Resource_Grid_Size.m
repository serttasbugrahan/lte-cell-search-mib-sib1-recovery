%buğrahan serttaş
%bu işlem ofdm demodülasyonu öncesi kaynak ızgarası uzunluğu için yapılır
function gridSize = lte_Resource_Grid_Size(enb)
normalCP = strcmpi(enb.CyclicPrefix, 'Normal');
gridSize = [12*enb.NDLRB, 14*normalCP + 12*(~normalCP),enb.CellRefP];%normal=72*14 extended = 72*12
end