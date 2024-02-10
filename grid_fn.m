%buğrahan serttaş
%boş kaynak ızgarası oluşturmak için kullanılır
function reGrid = grid_fn(enb)
normalCP = strcmpi(enb.CyclicPrefix, 'Normal');
d = [12*enb.NDLRB, 14*normalCP + 12*(~normalCP)];%normal=72*14 extended = 72*12
reGrid = zeros(d);
end