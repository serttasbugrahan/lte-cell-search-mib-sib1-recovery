%buğrahan serttaş
%PSS indis oluşturmak için yapılır PSS hesaplamalarından önce gerçekleşir
function PSS_ind = lte_PSS_Indices_fn(enb)
NscRB = 12;
normalCP = strcmpi(enb.CyclicPrefix, 'Normal');
L = 14*normalCP + 12*(~normalCP); % Number of symbols in one subframe
activeSubframeMod5 = 0;
l = (L/2) - 1;
if (mod(enb.NSubframe,5) == activeSubframeMod5)
    port =0;
    n = (0:61).';
    k = n - 31 + (NscRB*enb.NDLRB/2);
    k_len = numel(k);
    l = l*ones(k_len, 1);
    p = port.*ones(k_len, 1);
    PSS_ind = [k,l,p];
else
    PSS_ind = zeros(0,3);
end
% Convert 0based indices to 1based
PSS_ind = PSS_ind + 1;
% Extract k,l,p
k = PSS_ind(:,1);
l = PSS_ind(:,2);
p = PSS_ind(:,3);
opts.indexStyle='ind';
if (strcmpi(opts.indexStyle,'ind'))
    if (isempty(PSS_ind))
        PSS_ind = zeros(0,1);
    else
        normalCP = strcmpi(enb.CyclicPrefix, 'Normal');
        gridsize = [12*enb.NDLRB, 14*normalCP+12*(~normalCP)];
        % sub to ind conversion
        % Convert subscript based indices to linear indices
        PSS_ind = gridsize(1)*gridsize(2)*(p - 1) + gridsize(1)*(l - 1) + k;
        % gridsize(1) is the number of subcarriers, gridsize(2) is the
        % number of OFDM symbols in a subframe
    end
else
    % subs out, re-assemble k,l,p into matrix
    PSS_ind = [k,l,p];
end
end