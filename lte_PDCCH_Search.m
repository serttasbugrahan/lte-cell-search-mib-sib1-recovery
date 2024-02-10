%buğrahan serttaş
%PDCCH sonucunda bulunan DCI bitlerinin ayrıştırılması için yapılır
%PDCCH işleminden sonra uygulanır

function [decDCI,decDCIBits] = lte_PDCCH_Search(enb, pdcch, dciBits)
    % Validate and default the shared (mandatory and optional) parameters
    enbConfig = enb;
    ueConfig  = pdcch;
    % Input size validation
    if length(dciBits)<72
        error('lte:error','Input soft bits should be at least a CCE in length (72 bits).');
    end
    % Deduce total number of REGs (1 REG = 4 REs = 8 bits)
    enbConfig.NREG = floor(length(dciBits)/(72/9));
    
    % DCI formats for common and UE-specific search space
    %
    % Common search space DCI formats, note that Format1A is not listed
    % because it is the same size as Format0
    dciFormats{1} = {'Format0','Format1C'};
  
    % Cache the UE specific parameters and initially remove the RNTI since 
    % we will start by searching in the common search space
    pdcchConfig = ueConfig;
    pdcchConfig.ControlChannelType = 'PDCCH';  
    pdcchConfig = rmfield(pdcchConfig,'RNTI');
    pdcchConfig.SearchSpace = 'Common';
    % Ensure that 'DCIFormat' is not present in UE-specific parameters
    % since it may conflict with the use of lteDCI later 
    if isfield(pdcchConfig,'DCIFormat')
        pdcchConfig = rmfield(pdcchConfig,'DCIFormat');
    end   
    % Get DCI formats and lengths (relative to common search space 
    % for formats 0/1A and UE-specific otherwise), only the first format
    % for each unique message length is listed
    dciinfo = lte.internal.uniqueDCILengths(enbConfig,pdcchConfig);%NDBLR değerine göre dcı bit uzunluğu verilir
     
    % Identify UE-specific search space DCI formats
    dcimessages = fieldnames(dciinfo);
    uespecific_exclusionList = {'Format3','Format3A'};   % Remove 3/3A from this list
    dciFormats{2} = setdiff(dcimessages,uespecific_exclusionList);
    
    % PDCCH format for common search space can either be 2 or 3 (i.e.
    % aggregation level of 4 or 8 CCEs)
    startingPdcchFormat = 2;
       
    % Intermediate local variables
    idx = 1;
    decDCI = {};
    decDCIBits = {};
    reservedLoc = [];
   
    for searchType=1:2
        % UE-specific search space
        if(searchType == 2)
            pdcchConfig.RNTI = ueConfig.RNTI;
            pdcchConfig.SearchSpace = 'UESpecific';
            % Update the message lengths for UE-specific search space
            % (updates only relevant to formats 0/1A), only the first
            % format for each unique message length is listed
            dciinfo = lte.internal.uniqueDCILengths(enbConfig,pdcchConfig);
            % PDCCH format for UE-specific search space can be 0,1,2 or 3
            startingPdcchFormat = 0;
        end
        for pdcchFormat = 3:-1:startingPdcchFormat
            pdcchConfig.PDCCHFormat = pdcchFormat;

            % Performs common and/or ue-specific search depending on whether
            % the RNTI field is present in the channel/UE-specific parameters
            pdcchCandidates = ltePDCCHSpace(enbConfig,pdcchConfig,{'bits','1based'});%burayı anlayamadım

            % PDCCH candidates need not to be unique so picking the unique set
            % of candidates to optimize the search
            pdcchCandidates = unique(pdcchCandidates,'rows');

            pdcchCandidatesDims = size(pdcchCandidates);
            noOfCandidates = pdcchCandidatesDims(1);
            
            % Test each candidate in the search space for the presence of a decodable message format
            for candidate=1:noOfCandidates
                if(sum(reservedLoc == pdcchCandidates(candidate,1)/72) == 0)
                    if((pdcchCandidates(candidate,1)<length(dciBits)) && (pdcchCandidates(candidate,2)<=length(dciBits)))
                        input = dciBits(pdcchCandidates(candidate,1):pdcchCandidates(candidate,2));
                        for dciFormatIdx=1:length(dciFormats{searchType}) % Iterating through all DCI formats                     
                            df = dciFormats{searchType}{dciFormatIdx};        
                            mlength = dciinfo.(df);
                            [dciMessageBits,decRnti] = lteDCIDecode(mlength,input);
                            if(ueConfig.RNTI == decRnti && ~any(reservedLoc == pdcchCandidates(candidate,1)/72))
                                % Creating DCI message for successfully decoded PDCCH payload bits
                                dciMessage = struct;
                                dciMessageBits = transpose(dciMessageBits);
                                dciMessageBits = char('0' + dciMessageBits);
                                if isequal(dciMessageBits(1,1:1), '1')
                                    dciMessage.DCIFormat = 'Format1A';
                                elseif isequal(dciMessageBits(1,1:1), '0')
                                    dciMessage.DCIFormat;
                                end
                                dciMessage.CIF = bin2dec(num2str(dciMessageBits(2:2)));
                                dciMessage.AllocationType = bin2dec(num2str(dciMessageBits(3:3)));
                                dciMessage.Allocation = struct;
                                dciMessage.Allocation.Allocation = bin2dec(num2str(dciMessageBits(4:13)));
                                dciMessage.ModCoding = bin2dec(num2str(dciMessageBits(14:18)));
                                dciMessage.HARQNo = bin2dec(num2str(dciMessageBits(19:21)));
                                dciMessage.NewData = bin2dec(num2str(dciMessageBits(22:22)));
                                dciMessage.RV = bin2dec(num2str(dciMessageBits(23:24)));
                                dciMessage.TPCPUCCH = bin2dec(num2str(dciMessageBits(25:26)));
                                dciMessage.TDDIndex = bin2dec(num2str(dciMessageBits(27:27)));
                                dciMessage.SRSRequest = 0;
                                dciMessage.HARQACKResOffset = 0;
                                decDCIBits{idx} = dciMessageBits;
                                decDCI{idx} = dciMessage;
                                reservedLoc = [reservedLoc,(pdcchCandidates(candidate,1)/72:(pdcchCandidates(candidate,1)/72)+(2^pdcchFormat)-1)]; %#ok<AGROW>
                                idx = idx + 1;
                                break;     % Now break out of the format list loop since a format was found in the candidate
                            end
                        end
                    end
                end
            end
        end
    end
end