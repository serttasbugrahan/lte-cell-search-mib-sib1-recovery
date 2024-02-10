%buğrahan serttaş
%SSS indis değerlerini oluştumak için yapılır
%SSS işleminin başlangıcında kullanılır
function SSS_ind = lte_SSS_Indices_fn(enb)
NscRB = 12;
normalCP = strcmpi(enb.CyclicPrefix, 'Normal');
L = 14*normalCP + 12*(~normalCP); % Number of symbols in one subframe
activeSubframeMod5 = 0;
l = (L/2) - 2;
if (mod(enb.NSubframe,5) == activeSubframeMod5)
    port =0;
    n = (0:61).';
    k = n - 31 + (NscRB*enb.NDLRB/2);
    k_len = numel(k);
    l = l*ones(k_len, 1);
    p = port.*ones(k_len, 1);
    SSS_ind = [k,l,p];
else
    SSS_ind = zeros(0,3);
end
% Convert 0based indices to 1based
SSS_ind = SSS_ind + 1;
% Extract k,l,p
k = SSS_ind(:,1);
l = SSS_ind(:,2);
p = SSS_ind(:,3);
opts.indexStyle='ind';
if (strcmpi(opts.indexStyle,'ind'))
    if (isempty(SSS_ind))
        SSS_ind = zeros(0,1);
    else
        normalCP = strcmpi(enb.CyclicPrefix, 'Normal');
        gridsize = [12*enb.NDLRB, 14*normalCP+12*(~normalCP)];
        % sub to ind conversion
        % Convert subscript based indices to linear indices
        SSS_ind = gridsize(1)*gridsize(2)*(p - 1) + gridsize(1)*(l - 1) + k;
        % gridsize(1) is the number of subcarriers, gridsize(2) is the
        % number of OFDM symbols in a subframe
    end
else
    % subs out, re-assemble k,l,p into matrix
    SSS_ind = [k,l,p];
end
end