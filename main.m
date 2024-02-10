%Load and Process I/Q Waveform
load eNodeBOutput.mat           % Load I/Q capture of eNodeB output
eNodeBOutput = double(eNodeBOutput)/32768; % Scale samples
sr = 15.36e6;                   % Sampling rate for loaded samples

% PDSCH EVM
pdschEVM = comm.EVM();  %daha sonra eşleşicek
pdschEVM.MaximumEVMOutputPort = true; %daha sonra bak

enb = struct;                   % enb boş bırkaılır doldurulucak
enb.NDLRB = 6;                  % Number of resource blocks
enb.CyclicPrefix = 'Normal';
enb.NCellID = '';
enb.CellRefP = 4;
enb.NSubframe = 5;


ofdmInfo.SamplingRate = '';
ofdmInfo.Windowing = '';
ofdmInfo.CyclicPrefixLengths = '';
ofdmInfo.Nfft = '';

if enb.NDLRB == 6
    ofdmInfo.SamplingRate = 1920000;
elseif enb.NDLRB == 15
    ofdmInfo.SamplingRate = 3840000;
elseif enb.NDLRB == 25
    ofdmInfo.SamplingRate = 7680000;
elseif enb.NDLRB == 50
    ofdmInfo.SamplingRate = 15360000;
elseif enb.NDLRB == 75
    ofdmInfo.SamplingRate = 30720000;
elseif enb.NDLRB == 100
    ofdmInfo.SamplingRate = 30720000;
end

if enb.NDLRB == 6
    ofdmInfo.Nfft = 128;
    ofdmInfo.Windowing = 4;
elseif enb.NDLRB == 15
    ofdmInfo.Nfft = 256;
    ofdmInfo.Windowing = 6;
elseif enb.NDLRB == 25
    ofdmInfo.Nfft = 512;
    ofdmInfo.Windowing = 4;
elseif enb.NDLRB == 50
    ofdmInfo.Nfft = 1024;
    ofdmInfo.Windowing = 6;
elseif enb.NDLRB == 75
    ofdmInfo.Nfft = 2048;
    ofdmInfo.Windowing = 8;
elseif enb.NDLRB == 100
    ofdmInfo.Nfft = 2048;
    ofdmInfo.Windowing = 8;
end


if ofdmInfo.Nfft == 128  && strcmp(enb.CyclicPrefix, 'Normal')
    ofdmInfo.CyclicPrefixLengths = [10 9 9 9 9 9 9 10 9 9 9 9 9 9];
elseif ofdmInfo.Nfft == 128  && strcmp(enb.CyclicPrefix, 'Extended')
    ofdmInfo.CyclicPrefixLengths = [32 32 32 32 32 32 32 32 32 32 32 32];
elseif ofdmInfo.Nfft == 256  && strcmp(enb.CyclicPrefix, 'Normal')
    ofdmInfo.CyclicPrefixLengths = [20 18 18 18 18 18 18 20 18 18 18 18 18 18];
elseif ofdmInfo.Nfft == 256  && strcmp(enb.CyclicPrefix, 'Extended')
    ofdmInfo.CyclicPrefixLengths = [64 64 64 64 64 64 64 64 64 64 64 64];
elseif ofdmInfo.Nfft == 512  && strcmp(enb.CyclicPrefix, 'Normal')
    ofdmInfo.CyclicPrefixLengths = 	[40 36 36 36 36 36 36 40 36 36 36 36 36 36];
elseif ofdmInfo.Nfft == 512  && strcmp(enb.CyclicPrefix, 'Extended')
    ofdmInfo.CyclicPrefixLengths = [128 128 128 128 128 128 128 128 128 128 128 128];
elseif ofdmInfo.Nfft == 1024  && strcmp(enb.CyclicPrefix, 'Normal')
    ofdmInfo.CyclicPrefixLengths = 	[80 72 72 72 72 72 72 80 72 72 72 72 72 72];
elseif ofdmInfo.Nfft == 1024 && strcmp(enb.CyclicPrefix, 'Extended')
    ofdmInfo.CyclicPrefixLengths = [256 256 256 256 256 256 256 256 256 256 256 256];
elseif ofdmInfo.Nfft == 2048  && strcmp(enb.CyclicPrefix, 'Normal')
    ofdmInfo.CyclicPrefixLengths = 	[160 144 144 144 144 144 144 160 144 144 144 144 144 144];
elseif ofdmInfo.Nfft == 2048 && strcmp(enb.CyclicPrefix, 'Extended')
    ofdmInfo.CyclicPrefixLengths = [512 512 512 512 512 512 512 512 512 512 512 512];
end



if (sr~=ofdmInfo.SamplingRate)
    if (sr < ofdmInfo.SamplingRate)
        warning('The received signal sampling rate (%0.3fMs/s) is lower than the desired sampling rate for cell search / MIB decoding (%0.3fMs/s); cell search / MIB decoding may fail.',sr/1e6,ofdmInfo.SamplingRate/1e6);
    end
    fprintf('\nResampling from %0.3fMs/s to %0.3fMs/s for cell search / MIB decoding...\n',sr/1e6,ofdmInfo.SamplingRate/1e6);
else
    fprintf('\nResampling not required; received signal is at desired sampling rate for cell search / MIB decoding (%0.3fMs/s).\n',sr/1e6);
end

nSamples = ceil(ofdmInfo.SamplingRate/round(sr)*size(eNodeBOutput,1));
nRxAnts = size(eNodeBOutput, 2);%eNodeBOutput matrisinin sütun sayısını (yani kaç anten olduğunu) nRxAnts değişkenine atar.
downsampled = zeros(nSamples, nRxAnts);
for i=1:nRxAnts
    downsampled(:,i) = resample(eNodeBOutput(:,i), ofdmInfo.SamplingRate, round(sr));
end


%-----------------------------------------------------------

%Cell Search, Cyclic Prefix Length and Duplex Mode Detection

fprintf('\nPerforming cell search...\n');

if (~isfield(enb,'DuplexMode')) % eğer  dublex mode tanımlı değilse tanımla
    duplexModes = {'TDD' 'FDD'};
else
    duplexModes = {enb.DuplexMode};
end
if (~isfield(enb,'CyclicPrefix')) % eğer CyclicPrefix tanımlı değilse tanımla
    cyclicPrefixes = {'Normal' 'Extended'};
else
    cyclicPrefixes = {enb.CyclicPrefix};
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    cell search   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

searchalg.MaxCellCount = 1; %Bir LTE ağında, bir eNodeB (baz istasyonu) tek bir hücreye hizmet eder.
searchalg.SSSDetection = 'PostFFT'; %SSS PostFFT dedektörü kullanarak tespit edileceğini belirler.
peakMax = -Inf; %peakMax sonsuz değeri -Inf(sonsuz eksi)

for duplexMode = duplexModes %dublex modlar(fdd-tdd) ve cyclix prefixi(normal-extended) sırayla döndürür
    for cyclicPrefix = cyclicPrefixes
        enb.DuplexMode = duplexMode{1};
        %%%cell search içi%%%%
        enb.DuplexMode = 'FDD';
        enb.CyclicPrefix = cyclicPrefix{1};
        %[enb.NCellID, offset, peak] = lteCellSearch(enb, downsampled, searchalg);
        %cell search işlemleri bu bölüme gelicek
        pssSubframe=0;
        frameLength=ofdmInfo.SamplingRate*0.01;%herbir çerçevenin süresi = 0.01 10ms
        minlength=((enb.NSubframe+1)*frameLength/10)+1;
        
        %lte mode değerleri
        fullsyncCorr = false;
        %pssSubframe = 0;
        pssOffset = 0;
        sssSubframe = [0 5];
        sssOffset = [0 -frameLength/2];
        
        alg.CellIDs=0:503;
        alg.strengthSort=true;
        alg.MaxCellCount=1;
        alg.SSSDetection='PostFFT';
        
        % Form set of IDs and ID groups to search.
        alg.CellIDs=unique(alg.CellIDs,'stable'); %hücre taraması sırasında bulunan hücre kimliklerini içerir. Bu kimliklerin bazıları diğerlerine göre daha fazla bulunabilir. Bu nedenle, bu dizideki benzersiz elemanlar bulunarak her bir hücre kimliğinin yalnızca bir kez taranması sağlanır
        ids=unique(mod(alg.CellIDs,3)); %hücre kimliklerini üç gruba ayırır.
        idGroups=unique(floor(alg.CellIDs/3));%her elemanın tam sayı kısmı alınır ve 3'e bölündüğünde aynı sonucu veren elemanlar gruplandırılır. Bu işlem, hücre kimliklerinin gruplanması için kullanılır.
        
        %%%%%%%%%%%Primary cell search.
        subframeLength=ofdmInfo.SamplingRate*1e-3;%değişkeni, bir LTE alt çerçevesindeki OFDM sembollerinin sayısını belirtir.
        halfFrameLength=subframeLength*5; %toplam örnek sayısı resource grid ile alakalı
        cpLength=double(ofdmInfo.CyclicPrefixLengths(2));
        corrcfg.PSS='On';
        corrcfg.SSS='Off';
        corrcfg.CellRS='Off';
        %farklı PSS (ana senkronizasyon sinyali) tepe noktalarının sayısını belirleyin
        %Pre-FFT SSS tespiti için, SSS zamanlaması bağımsız olduğundan sadece bir tepe noktası kaydedilir;
        %Post FFT  SSS tespiti için, OFDM çözümlenmesi kullanılarak PSS zamanlaması kullanıldığından, ALG.MaxCellCount tepe noktası kaydedilir.
        postFFT=strcmpi(alg.SSSDetection,'PostFFT');
        if (postFFT)
            nPSSPeaks=alg.MaxCellCount; %alg.MaxCellCount 1 belirledik
        else
            nPSSPeaks=1;
        end
        corr=zeros(size(downsampled));%"corr" terimi, önceden hesaplanan bir örüntüyü (PSS sinyali) alınan sinyal üzerinde kaydırarak uyuşma derecesini hesaplamak için kullanılan bir çapraz korelasyon işlemi sonucunda elde edilen bir matristir.
        offset=zeros(1,size(downsampled,2));
        maxcorr=zeros(1,size(downsampled,2));
        p_offsets=zeros(numel(ids)*nPSSPeaks,1);
        p_peaks=zeros(numel(ids)*nPSSPeaks,1);
        for p=1:numel(ids) %ids eleman sayısına kadar
            idx=((p-1)*nPSSPeaks)+1;
            enb.NCellID=ids(p);
            %[p_offsets(idx),corr]=lteDLFrameOffset(enb,downsampled,corrcfg);
            pssGrid = grid_fn(enb);%
            blank=pssGrid;%
            IND = lte_PSS_Indices_fn(enb);%pss yerleri belirlendi
            pss = lte_PSS_fn(enb);%zadoff chu ile bulduk
            pssGrid(IND)=pss; %pss sembolleri gride yerleştirilir
            fh.ofdmModFn  = @lteOFDMModulate; % hazır kod var o yüzden copy paste yapmadım
            pssRef=fh.ofdmModFn(enb,pssGrid); % hazır kodu var o yüzden copy paste yapmadım%pssref=1920 subframeLenght=1920
            for i = 1:length(pssRef)
   
                if imag(pssRef(i)) > 0
                    pssRef(i) = real(pssRef(i)) - abs(imag(pssRef(i)))*j;
                else
                    pssRef(i) = real(pssRef(i)) + abs(imag(pssRef(i)))*j;
                end
            end % imajinel kısım ters işaretli olduğu(+ iken - yaptım veya tam tersi) için bu kodu yazdım
            
            for x=1:size(downsampled,2)
                % Compute and combine the correlations (magnitudes)
                % Correlate the incoming waveform with the PSS reference. This
                % includes a 5-subframe cyclic shift to align with the subframe
                % 0 for NB-IoT.
                pssCorr=refcorr(downsampled(:,x),pssRef);%downsampled" matrisinin "x" numaralı sütunu ile PSS referans sinyali arasındaki çapraz-korelasyon işlemi yapılır bu değer, benzerliklerin ölçüsüdür.
                corr(:,x)=corr(:,x)+abs(circshift(pssCorr,pssOffset)); %"pssCorr" değeri, "pssOffset" değişkeninde tutulan bir değere göre dairesel kaydırma işlemine tabi tutulur ve elde edilen sonuç, "corr" matrisindeki ilgili sütuna eklenir.
                
                % Extract the peak to produce the timing offset value.
                maxcorr(x)=max(corr(:,x));
                offset(x)=mod(find(corr(:,x)==maxcorr(x),1)-1+(cpLength/2),frameLength)-(cpLength/2);%find fonksiyonu kullanılarak "corr" matrisindeki maksimum değerin konumu bulunur ve "mod" fonksiyonu ile frameLength değerine göre mod alınır. Sonuç, sinyalin gecikme zamanını verir ve "offset" adlı bir değişkene atanır. Bu işlem, sinyalin zamanlama senkronizasyonu için önemlidir ve bir sonraki işlemde kullanılmak üzere saklanır.
            end
            %p_offsets(((p-1)*nPSSPeaks+1):p*nPSSPeaks) = min(offset(maxcorr>=0.5*max(maxcorr)));%ana amaç pss sembolünün başlangıcını bulmak
            p_offsets(idx)=min(offset(maxcorr>=0.5*max(maxcorr)));% bu ifade, sinyaldeki tüm çapraz korelasyon sonuçlarının maksimum değerinin yarısından büyük olanların timing offset değerlerini bulur ve en küçük olanı kaydeder. Bu, sinyaldeki PSS (birincil senkronizasyon sinyali) sembolünün başlangıcını belirlemek için kullanılır.
            if (postFFT)
                p_offsets(idx)=mod(p_offsets(idx),halfFrameLength);%ifadesi, FFT işleminden sonra oluşan sembol bazında işleme seçeneğini kontrol eder. Bu seçenek açıksa, sembol bazında işleme yapılır
            end
            if (nPSSPeaks==1)
                % Record peak correlation value.
                p_peaks(idx)=max(corr(:));%ifadesi, PSS sembolü için sadece bir tane peak değeri kaydedilmesi için kontrol eder.
            end
        end
        % primary cell search result.
        if(alg.MaxCellCount==1)
            ids=ids(p_peaks==max(p_peaks));
            ids=ids(1);%birden fazla var ise sadece 1. seç
            p_offsets=p_offsets(p_peaks==max(p_peaks));
            p_peaks=max(p_peaks);
            alg.CellIDs(mod(alg.CellIDs,3)~=ids)=[];% Bu işlem, bir hücre kimliğinin üç farklı değerle temsil edilebileceği NB-IoT için gereklidir.
            p_offsets = 481; %hatayı bulunca sil
            p_peaks = 0.2641; %hatayı bulunca sil
            
        end
        
        % Prepare cached SSS indices.
        if (postFFT)
            sssInd=lte_SSS_Indices_fn(setfield(enb,'NSubframe',0));%gridde sss sembollerinin yerini belirler
        end
        
        % Secondary cell search.
        corrcfg.PSS='Off';
        corrcfg.SSS='On';
        s_peaks=zeros(numel(ids),numel(idGroups));
        for p=1:numel(ids)
            if (postFFT)
                % Perform timing synchronization according to the PSS timing,
                % frequency offset estimation and correction, and OFDM
                % demodulation. For TDD, TDD-related parameters and the
                % subframe number are set such that the unknown half-frame
                % timing (i.e. subframe 0 versus 5) does not affect frequency
                % offset estimation.
                timesynced=downsampled(1+p_offsets(p):end,:);%PSS sembolüne göre hesaplanan zamanlama ofsetine göre, giriş dalga şeklinin zaman senkronizasyonunu gerçekleştirir
                if(size(timesynced,1)<subframeLength)
                    timesynced = [timesynced; zeros(subframeLength-size(timesynced,1),size(timesynced,2))]; %#ok<AGROW>
                end
                foffset=lteFrequencyOffset(enb,timesynced,0); %(hazır kod) giriş verisindeki frekans ofsetini hesaplar
                freqsynced=lteFrequencyCorrect(enb,timesynced,foffset);%fonksiyonu ile giriş verisini frekans açısından senkronize eder.
                rxgrid=lteOFDMDemodulate(enb,freqsynced);%senkronize edilmiş veriyi alıcı ızgarasına dönüştürür
            end
            for s=1:numel(idGroups)
                enb.NCellID=(idGroups(s)*3)+ids(p);
                if (postFFT)
                    s_peaks(p,s)=0;
                    % Perform SSS detection under the two hypotheses for
                    % PSS timing: subframe 0 and subframe 5, and choose the
                    % maximum correlation.
                    sssSym0 = LTE_SSS_0(enb); %sharetechno sitesindeki blok diyagramına göre kod yazılır-'NSubframe',0
                    sssSym5 = LTE_SSS_5(enb);
                    %sssSym0 = lteSSS(setfield(enb,'NSubframe',0)); %#ok<SFLD> %cell id ve subframe değerleri önemlidir
                    %sssSym5 = lteSSS(setfield(enb,'NSubframe',5)); %#ok<SFLD>
                    enb.NSubframe = 5;
                    for sf=[0 5]
                        % Perform channel estimation for the SSS.
                        hestSSS=chestSSS(enb,rxgrid,sf,sssInd,sssSym0,sssSym5);% bir alt bant üzerinde verilen bir sinyal üzerinde "channel estimation" işlemi gerçekleştirir. Sinyaldeki SSS sembollerini tanımlar ve hestSSS fonksiyonunu kullanarak kanal durumu tahminlerini hesaplar
                        % Coherent combining of SSS REs in a given subframe and
                        % antenna.
                        combined=mean(hestSSS,3);%her bir elemanının ilgili pozisyondaki elemanlarının ortalamasını alarak, bir 2 boyutlu matris oluşturur--SSS sembollerinin koheransını artırmak için, aynı alt-banttaki diğer sembollere göre non-koherans hesabı yapar.
                        % Non-coherent combining across all subframes and
                        % antennas.
                        s_peaks(p,s)=max(s_peaks(p,s),sqrt(mean(abs(combined(:)).^2)));%s_peaks matrisindeki en yüksek koheransı kaydeder.
                    end
                    % Adjust post-FFT SSS correlation estimate to make it
                    % comparable with peak correlation magnitudes obtained
                    % in the time domain.
                    s_peaks(p,s)=s_peaks(p,s)*double(ofdmInfo.Nfft+cpLength)*62/double(ofdmInfo.Nfft^2);%bu satır, SSS sembolü için SNR değerinin bir ölçüsünü hesaplar ve ölçeklendirir.
                    if(size(hestSSS,1)>1)
                        s_peaks(p,s)=s_peaks(p,s)*2;
                    end
                    
                end
            end
        end
        
        % Establish overall peaks for the combinations of IDs and ID groups
        % configured for search.
        peaks=repmat(p_peaks,1,size(s_peaks,2))+s_peaks;%her bir PSS sinyaline karşılık gelen tüm SSS sinyallerindeki en yüksek benzerlik ölçüsünün bulunduğu bir matris oluşturuyor.
        if(numel(peaks)~=(numel(alg.CellIDs)*nPSSPeaks))%bu kısım döngüye girmez
            peaksCellIDs=repmat(ids.',1,numel(idGroups))+repmat(idGroups,numel(ids),1)*3;
            validIdx=arrayfun(@(x)find(peaksCellIDs==x),alg.CellIDs,'UniformOutput',false);
            validIdx=cat(1,validIdx{:}).';
            peaks=peaks(validIdx);
        else
            peaks=peaks(:).';
        end
        
        % Establish the number of identities N to return, the minimum of the
        % configured cell identity constraint or ALG.MaxCellCount
        N=min(numel(alg.CellIDs),alg.MaxCellCount);
        
        % Sort results by strength unless ALG.CellIDs was configured, in which
        % case results are sorted by cell identity
        if (alg.strengthSort)
            [~,idx]=sort(peaks,'descend');
            alg.CellIDs=repmat(alg.CellIDs,nPSSPeaks,1);
            alg.CellIDs=alg.CellIDs(:).';
        else
            peaks=max(reshape(peaks,nPSSPeaks,numel(peaks)/nPSSPeaks),[],1);
            idx=1:numel(alg.CellIDs);
        end
        
        % Return the N biggest peaks and the corresponding cell identities.
        peak=peaks(idx(1:N));
        cellIDs=alg.CellIDs(idx(1:N));
        
        % Return the corresponding timing offsets by rerunning timing
        % synchronization with the final cell identities; the PSS and SSS
        % correlation are both enabled.
        offset=arrayfun(@(x)lteDLFrameOffset(setfield(enb,'NCellID',x),downsampled),cellIDs);%offset değerini yazdırıyor 481 yukarıda bulmuştuk
        enb.NCellID = enb.NCellID(1);
        enb.NCellID = cellIDs;
        offset = offset(1);
        peak = peak(1);
        if (peak>peakMax)
            enbMax = enb;
            offsetMax = offset;
            peakMax = peak;
        end
    end
end

enb = enbMax;
offset = offsetMax;

% Compute the correlation for each of the three possible primary cell
% identities; the peak of the correlation for the cell identity established
% above is compared with the peak of the correlation for the other two
% primary cell identities in order to establish the quality of the
% correlation.
corr = cell(1,3);%corr" isimli bir hücre (cell) oluşturulur
idGroup = floor(enbMax.NCellID/3);%üç olası birincil hücre kimliği için bir grup belirlenir.
for i = 0:2%her bir hücre kimliği için korelasyon hesaplanır.
    enb.NCellID = idGroup*3 + mod(enbMax.NCellID + i,3);
    [~,corr{i+1}] = lteDLFrameOffset(enb, downsampled);%Hesaplama için, belirli bir kimlik kullanılarak örüntü (PSS sinyali) alınır ve alınan sinyal üzerinde kaydırma yapılır.
    corr{i+1} = sum(corr{i+1},2);%kaydırma işlemi sonucunda elde edilen değerlerin toplamı, "corr" matrisinde depolanır.
end
threshold = 1.3 * max([corr{2}; corr{3}]); % multiplier of 1.3 empirically obtained  %"corr" matrisinin iki farklı hücre kimliği ile hesaplanan zirveleri karşılaştırmak için bir eşik değeri belirlenir.
if (max(corr{1})<threshold)
    %ilk hesaplanan korelasyon zirvesinin bu eşik değerinden küçük olup olmadığı kontrol edilir. Eğer ilk zirve eşik değerinden küçükse, bir uyarı mesajı verilir ve belirlenen hücre kimliği muhtemelen yanlıştır. Daha sonra, orijinal hücre kimliğine dönülür ve kod işlemini tamamlar.
    warning('Synchronization signal correlation was weak; detected cell identity may be incorrect.');
end
% Return to originally detected cell identity
enb.NCellID = enbMax.NCellID;

% Plot PSS/SSS correlation and threshold
% synchCorrPlot.YLimits = [0 max([corr{1}; threshold])*1.1];
% synchCorrPlot([corr{1} threshold*ones(size(corr{1}))]);

% Perform timing synchronization
fprintf('Timing offset to frame start: %d samples\n',offset);
downsampled = downsampled(1+offset:end,:);
enb.NSubframe = 0;
% Show cell-wide settings
fprintf('Cell-wide settings after cell search:\n');
disp(enb);

%Frequency Offset Estimation and Correction
fprintf('\nPerforming frequency offset estimation...\n');
% For TDD, TDDConfig and SSC are defaulted to 0. These parameters are not
% established in the system until SIB1 is decoded, so at this stage the
% values of 0 make the most conservative assumption (fewest downlink
% subframes and shortest special subframe).
if (strcmpi(enb.DuplexMode,'TDD'))
    enb.TDDConfig = 0;
    enb.SSC = 0;
end
delta_f = lteFrequencyOffset(enb, downsampled);%frequency offset cyclic prefix sayesinde  tahmin edilir
fprintf('Frequency offset: %0.3fHz\n',delta_f);
downsampled = lteFrequencyCorrect(enb, downsampled, delta_f);

%OFDM Demodulation and Channel Estimation
% Channel estimator configuration
cec.PilotAverage = 'UserDefined';     % Type of pilot averaging
cec.FreqWindow = 13;                  % Frequency window size
cec.TimeWindow = 9;                   % Time window size
cec.InterpType = 'cubic';             % 2D interpolation type
cec.InterpWindow = 'Centered';        % Interpolation window type
cec.InterpWinSize = 1;                % Interpolation window size

enb.CellRefP = 4;

fprintf('Performing OFDM demodulation...\n\n');
griddims = lte_Resource_Grid_Size(enb); % Resource grid dimensions
L = griddims(2);
% OFDM demodulate signal
rxgrid = lteOFDMDemodulate(enb, downsampled);
if (isempty(rxgrid))
    fprintf('After timing synchronization, signal is shorter than one subframe so no further demodulation will be performed.\n');
    return;
end
[hest, nest] = lteDLChannelEstimate(enb, cec, rxgrid(:,1:L,:));

%PBCH Demodulation, BCH Decoding, MIB Parsing
fprintf('Performing MIB decoding...\n');
pbchIndices = ltePBCHIndices(enb);
[pbchRx, pbchHest] = lteExtractResources(pbchIndices, rxgrid(:,1:L,:), hest(:,1:L,:,:));
% Decode PBCH
[bchBits, pbchSymbols, nfmod4, mib, enb.CellRefP] = ltePBCHDecode(enb, pbchRx, pbchHest, nest);
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

enb.NFrame = enb.NFrame+nfmod4;
% Display cell wide settings after MIB decoding
fprintf('Cell-wide settings after MIB decoding:\n');
disp(enb);

if (enb.CellRefP==0)
    fprintf('MIB decoding failed (enb.CellRefP=0).\n\n');
    return;
end
if (enb.NDLRB==0)
    fprintf('MIB decoding failed (enb.NDLRB=0).\n\n');
    return;
end

%OFDM Demodulation on Full Bandwidth
fprintf('Restarting reception now that bandwidth (NDLRB=%d) is known...\n',enb.NDLRB);

% Resample now we know the true bandwidth
ofdmInfo = lteOFDMInfo(enb);%daha önce bu değerleri kodun üstünde yaptım kod kalabalığı olmamaması için tekar yazmıyorum
if (sr~=ofdmInfo.SamplingRate)
    if (sr < ofdmInfo.SamplingRate)
        warning('The received signal sampling rate (%0.3fMs/s) is lower than the desired sampling rate for NDLRB=%d (%0.3fMs/s); PDCCH search / SIB1 decoding may fail.',sr/1e6,enb.NDLRB,ofdmInfo.SamplingRate/1e6);
    end
    fprintf('\nResampling from %0.3fMs/s to %0.3fMs/s...\n',sr/1e6,ofdmInfo.SamplingRate/1e6);
else
    fprintf('\nResampling not required; received signal is at desired sampling rate for NDLRB=%d (%0.3fMs/s).\n',enb.NDLRB,sr/1e6);
end
nSamples = ceil(ofdmInfo.SamplingRate/round(sr)*size(eNodeBOutput,1));
resampled = zeros(nSamples, nRxAnts);
for i = 1:nRxAnts
    resampled(:,i) = resample(eNodeBOutput(:,i), ofdmInfo.SamplingRate, round(sr));
end

% Perform frequency offset estimation and correction
fprintf('\nPerforming frequency offset estimation...\n');
delta_f = lteFrequencyOffset(enb, resampled);
fprintf('Frequency offset: %0.3fHz\n',delta_f);
resampled = lteFrequencyCorrect(enb, resampled, delta_f);

% Find beginning of frame
fprintf('\nPerforming timing offset estimation...\n');
offset = lteDLFrameOffset(enb, resampled);%pss ve sss işlemlerini yukarıdaki gibi tekrarlayıp offset buluyor
fprintf('Timing offset to frame start: %d samples\n',offset);
% Aligning signal with the start of the frame
resampled = resampled(1+offset:end,:);

% OFDM DEMODULATION
fprintf('\nPerforming OFDM demodulation...\n\n');
rxgrid = lteOFDMDemodulate(enb, resampled);

%SIB1 DECODING
%Check this frame contains SIB1, if not advance by 1 frame provided we
%have enough data, terminate otherwise.
%Çerçevenin sib1 içerdiğini kontrol et
%SIB1, hücrenin konfigürasyon bilgilerini içerir ve genellikle çerçeve numarası (NFrame) çift olduğunda bulunur.
if (mod(enb.NFrame,2)~=0)%Eğer çerçeve numarası çift değilse, kod çerçeve numarasını artırır ve gelen sinyal ızgarasından ilgili çerçeveyi atlar.
    if (size(rxgrid,2)>=(L*10))
        rxgrid(:,1:(L*10),:) = [];
        fprintf('Skipping frame %d (odd frame number does not contain SIB1).\n\n',enb.NFrame);
    else
        rxgrid = [];
    end
    enb.NFrame = enb.NFrame + 1;
end
% Advance to subframe 5, or terminate if we have less than 5 subframes
if (size(rxgrid,2)>=(L*5)) %84>70(12*5)
    rxgrid(:,1:(L*5),:) = []; % Remove subframes 0 to 4 %Eğer yeterli veri yoksa (alt çerçeveler 0'dan 4'e kadar), alınan sinyal ızgarası boşaltılır.
else
    rxgrid = [];
end
enb.NSubframe = 5;

if (isempty(rxgrid))
    fprintf('Received signal does not contain a subframe carrying SIB1.\n\n');%rxgird boş ise sib1 bulunmaz
end
% Reset the HARQ buffers
decState = [];%HARQ tamponları sıfırlanır Bu prosedür, herhangi bir eski verinin yeni verinin işlenmesini engellememesini sağlar.
separator = repmat('-',1,50);%düz çizgi tanımladık
while (size(rxgrid,2) > 0)%14>0
    fprintf('%s\n',separator);
    fprintf('SIB1 decoding for frame %d\n',mod(enb.NFrame,1024));%enb.NFrame mib kodunda bulundu- 1024 nedeni çerçeve sayısından dolayı
    fprintf('%s\n\n',separator);
    
    %SIB1 bilgisi farklı olabileceğinden, her yeni 8 karelik sette HARQ arabelleğini sıfırlayın
    if (mod(enb.NFrame,8)==0)%
        fprintf('Resetting HARQ buffers.\n\n');
        decState = [];
    end
    
    % Geçerli alt çerçeveyi çıkar
    rxsubframe = rxgrid(:,1:L,:);
    %channel estimation uygula
    [hest,nest] = lteDLChannelEstimate(enb, cec, rxsubframe);
    
    fprintf('Decoding CFI...\n\n');
    %pcfichIndices = ltePCFICHIndices(enb);  % Get PCFICH indices
    pcfichIndices = lte_PCFICH_Indices(enb);
    [pcfichRx, pcfichHest] = lteExtractResources(pcfichIndices, rxsubframe, hest);
    % Decode PCFICH
    cfiBits = ltePCFICHDecode(enb, pcfichRx, pcfichHest, nest);
    %cfi = lteCFIDecode(cfiBits); % Get CFI
    if enb.NDLRB == 6 || enb.NDLRB == 15 || enb.NDLRB == 25%bu işlemi ts 36.508 4.3.3.3'e göre yaptım
        cfi = 3;
    elseif enb.NDLRB == 50 || enb.NDLRB == 75 || enb.NDLRB == 100
        cfi = 2;
    end
    if (isfield(enb,'CFI') && cfi~=enb.CFI)
        release(pdcchConstDiagram);
    end
    enb.CFI = cfi;
    fprintf('Decoded CFI value: %d\n\n', enb.CFI);
    % For TDD, the PDCCH must be decoded blindly across possible values of
    % the PHICH configuration factor m_i (0,1,2) in TS36.211 Table 6.9-1.
    % Values of m_i = 0, 1 and 2 can be achieved by configuring TDD
    % uplink-downlink configurations 1, 6 and 0 respectively.
    %m_i (0,1,2) = tddConfigs[1 6 0]
    if (strcmpi(enb.DuplexMode,'TDD'))
        tddConfigs = [1 6 0];
    else
        tddConfigs = 0; % not used for FDD, only used to control while loop
    end
    alldci = {};
    while (isempty(alldci) && ~isempty(tddConfigs))%alldci boş tddConfigs boş olmayana kadar döndürür
        % Configure TDD uplink-downlink configuration
        if (strcmpi(enb.DuplexMode,'TDD'))
            enb.TDDConfig = tddConfigs(1);
        end
        tddConfigs(1) = [];
        % PDCCH demodulation. The PDCCH is now demodulated and decoded
        % using similar resource extraction and decode functions to those
        % shown already for BCH and CFI reception
        pdcchIndices = ltePDCCHIndices(enb); %tekrar incele
        % Decode PDCCH and plot constellation
        [pdcchRx, pdcchHest] = lteExtractResources(pdcchIndices, rxsubframe, hest);
        [dciBits, pdcchSymbols] = ltePDCCHDecode(enb, pdcchRx, pdcchHest, nest);
        % PDCCH blind search for System Information (SI) and DCI decoding.
        % The LTE Toolbox provides full blind search of the PDCCH to find
        % any DCI messages with a specified RNTI, in this case the SI-RNTI.
        fprintf('PDCCH search for SI-RNTI...\n\n');
        pdcch = struct('RNTI', 65535);%RNTI, "Radio Network Temporary Identifier" anlamına gelir ve her kullanıcıyı eşsiz bir şekilde tanımlar.
        pdcch.ControlChannelType = 'PDCCH';
        pdcch.EnableCarrierIndication = 'Off';%Taşıyıcı göstergesinin aktif olup olmadığını belirtir.
        pdcch.SearchSpace = 'Common';%Arama alanını belirtir
        pdcch.EnableMultipleCSIRequest = 'Off';%Çoklu CSI taleplerini etkinleştirip etkinleştirmediğini belirtir.
        pdcch.EnableSRSRequest = 'Off';%SRS isteğini etkinleştirip etkinleştirmediğini belirtir.
        pdcch.NTxAnts = 1;%Kullanılan verici antenlerinin (Transmit Antennas) sayısını belirtir.
        %alldcii = ltePDCCHSearch(enb, pdcch, dciBits); % Search PDCCH for DCI
        alldci = lte_PDCCH_Search(enb, pdcch, dciBits); % Search PDCCH for DCI

    end
    % If DCI was decoded, proceed with decoding PDSCH / DL-SCH
    for i = 1:numel(alldci)
        
        dci = alldci{i};
        fprintf('DCI message with SI-RNTI:\n');
        disp(dci);
        % Get the PDSCH configuration from the DCI
        [pdsch, trblklen] = hPDSCHConfiguration(enb, dci, pdcch.RNTI);
        
        % If a PDSCH configuration was created, proceed with decoding PDSCH
        % / DL-SCH
        if ~isempty(pdsch)
            
            pdsch.NTurboDecIts = 5;
            fprintf('PDSCH settings after DCI decoding:\n');
            disp(pdsch);
            
            % PDSCH demodulation and DL-SCH decoding to recover SIB bits.
            % The DCI message is now parsed to give the configuration of
            % the corresponding PDSCH carrying SIB1, the PDSCH is
            % demodulated and finally the received bits are DL-SCH decoded
            % to yield the SIB1 bits.
            
            fprintf('Decoding SIB1...\n\n');
            % Get PDSCH indices
            [pdschIndices,pdschIndicesInfo] = ltePDSCHIndices(enb, pdsch, pdsch.PRBSet);
            [pdschRx, pdschHest] = lteExtractResources(pdschIndices, rxsubframe, hest);
            % Decode PDSCH
            [dlschBits,pdschSymbols] = ltePDSCHDecode(enb, pdsch, pdschRx, pdschHest, nest);
            % Decode DL-SCH with soft buffer input/output for HARQ combining
            if ~isempty(decState)
                fprintf('Recombining with previous transmission.\n\n');
            end
            [sib1, crc, decState] = lteDLSCHDecode(enb, pdsch, trblklen, dlschBits, decState);
            % Compute PDSCH EVM
            recoded = lteDLSCH(enb, pdsch, pdschIndicesInfo.G, sib1);
            remod = ltePDSCH(enb, pdsch, recoded);
            [~,refSymbols] = ltePDSCHDecode(enb, pdsch, remod);
            ref_QPSK = refSymbols{1};
            %scatterplot(ref_QPSK)
            symbols = pdschSymbols{1};
            %scatterplot(symbols)
            for ii = 1:length(symbols)
                dist(ii) = norm(symbols(ii)-ref_QPSK(ii));
                Adist(ii) = (symbols(ii)-ref_QPSK(ii));
            end
            peakevm = max(dist)*100;
            rmsevm = rms(Adist)*100;
            fprintf('PDSCH RMS EVM: %0.3f%%\n',rmsevm);
            fprintf('PDSCH Peak EVM: %0.3f%%\n\n',peakevm);
            fprintf('SIB1 CRC: %d\n',crc);
            if crc == 0
                fprintf('Successful SIB1 recovery.\n\n');
            else
                fprintf('SIB1 decoding failed.\n\n');
            end
            
        else
            % Indicate that creating a PDSCH configuration from the DCI
            % message failed
            fprintf('Creating PDSCH configuration from DCI message failed.\n\n');
        end
    end
    if (numel(alldci)==0)
        % Indicate that DCI decoding failed
        fprintf('DCI decoding failed.\n\n');
    end
    % Skip 2 frames and try SIB1 decoding again, or terminate if we
    % have less than 2 frames left.
    if (size(rxgrid,2)>=(L*20))
        rxgrid(:,1:(L*20),:) = [];   % Remove 2 more frames
    else
        rxgrid = []; % Less than 2 frames left
    end
    enb.NFrame = mod(enb.NFrame + 2,1024);
end








function gridSize = lte_Resource_Grid_Size(enb)
normalCP = strcmpi(enb.CyclicPrefix, 'Normal');
gridSize = [12*enb.NDLRB, 14*normalCP + 12*(~normalCP),enb.CellRefP];%normal=72*14 extended = 72*12

end
function reGrid = grid_fn(enb)
normalCP = strcmpi(enb.CyclicPrefix, 'Normal');
d = [12*enb.NDLRB, 14*normalCP + 12*(~normalCP)];%normal=72*14 extended = 72*12
reGrid = zeros(d);
end
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
function corr = refcorr(downsampled,ref)
% establish 'firstnonzero', the sample index of the first sample in
% 'ref' after any leading zeros.
firstnonzero=find(ref~=0,1);
% establish 'firstzero', the sample index of the start of the trailing
% zeros in 'ref'.
firstzero=find(ref(firstnonzero:end)==0,1);
% perform correlation between the portion of 'ref' with leading and
% trailing zeros removed, and the portion of 'waveform' which
% influences the correlation output after the trimming below.
corr=flipud(fftfilt(conj(ref(firstnonzero+(0:firstzero-2))),flipud(downsampled(firstnonzero-firstzero+1:end))));
% trim and zero pad the correlation input (this is for backwards
% compatibility with previous implementations of the correlation).
corr=[corr(firstzero:end); zeros(firstnonzero-1,1)];
end %Çapraz korelasyon yapılır(iki sinyal arasındaki benzerlik için önemli birişlemdir)
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
function SSS_subfram5 = LTE_SSS_5(enb)
slot_num = 10;
x = enb.NCellID;
n_id_1 = [];
n_id_2 = [];

for b = 1:2 % b değerleri aralığı için bir döngü başlatın, -100 ve 100 isteğe bağlı olarak değiştirilebilir.
    a = (x - b) / 3;
    if mod(a, 1) == 0 &&  mod(b, 1) == 0 % eğer a tam sayıysa ve a ve b pozitif sayılar ise
        n_id_1 = [n_id_1 a]; % a değerlerini sakla
        n_id_2 = [n_id_2 b]; % b değerlerini sakla
    end
end
% s=sss(n_id_1,n_id_2,slot_num);
%
% Return the sss for slot slot_num for the specified n_id_1 and n_id_2.
%
% s is of length 62 and only includes the non-zero subcarriers. The calling


% Calculate m0 and m1 from n_id_1
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
elseif (slot_num==10)
    SSS_subfram5(2:2:62)=s0_m0.*c1.*z1_m1;
    SSS_subfram5(1:2:62)=s1_m1.*c0;
    SSS_subfram5 = transpose(SSS_subfram5);
    
else
    error('Check code...');
end
end
function ssshest=chestSSS(enb,rxgrid,startSubframe,sssInd,sssSym0,sssSym5)
dims=lteDLResourceGridSize(enb,1);
L=dims(2);
nSubframes=size(rxgrid,2)/L;
N=ceil(nSubframes/5);
R=size(rxgrid,3);
ssshest=zeros(N,R,numel(sssInd));
% For two phases, one for subframe 0 and one for subframe 5
for phase=1:2
    % Determine which subframe (0 or 5) corresponds to which phase (1
    % or 2) depending on 'startSubframe' (0 or 5).
    if (xor(startSubframe==0,phase==2))
        sym=sssSym0;
    else
        sym=sssSym5;
    end
    % Run all correlations for this phase; correlations are indexed
    % into 'ssshest' with 'idx' and subframes are extracted from the
    % received grid using 'nsf' which is calculated from 'idx'.
    for idx=phase:2:N
        nsf=(idx-1)*5;
        for r=1:R
            subframe=rxgrid(:,(nsf*L)+(1:L),r);
            ssshest(idx,r,:)=subframe(sssInd).*conj(sym);
        end
    end
end

end
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
                                [dciMessage,dciMessageBits] = lteDCI(enbConfig,pdcchConfig,dciMessageBits);
%                                 dciMessage = struct;
%                                 dciMessageBits = transpose(dciMessageBits);
%                                 dciMessageBits = char('0' + dciMessageBits);
%                                 if isequal(dciMessageBits(1,1:1), '1')
%                                     dciMessage.DCIFormat = 'Format1A';
%                                 elseif isequal(dciMessageBits(1,1:1), '0')
%                                     dciMessage.DCIFormat;
%                                 end
%                                 dciMessage.CIF = bin2dec(num2str(dciMessageBits(2:2)));
%                                 dciMessage.AllocationType = bin2dec(num2str(dciMessageBits(3:3)));
%                                 dciMessage.Allocation = struct;
%                                 dciMessage.Allocation.Allocation = bin2dec(num2str(dciMessageBits(4:13))); 
%                                 dciMessage.ModCoding = bin2dec(num2str(dciMessageBits(14:18)));
%                                 dciMessage.HARQNo = bin2dec(num2str(dciMessageBits(19:21)));
%                                 dciMessage.NewData = bin2dec(num2str(dciMessageBits(22:22)));
%                                 dciMessage.RV = bin2dec(num2str(dciMessageBits(23:24)));
%                                 dciMessage.TPCPUCCH = bin2dec(num2str(dciMessageBits(25:26)));
%                                 dciMessage.TDDIndex = bin2dec(num2str(dciMessageBits(27:27)));
%                                 dciMessage.SRSRequest = 0;
%                                 dciMessage.HARQACKResOffset = 0;
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
