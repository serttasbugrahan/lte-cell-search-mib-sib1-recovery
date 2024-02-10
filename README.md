# Cell Search, MIB and SIB1 Recovery

This example shows how to fully synchronize, demodulate and decode a live eNodeB signal by using Matlab. Before the user equipment (UE) can communicate with the network, it must perform cell search and selection procedures, and then obtain initial system information. This process involves acquiring slot and frame synchronization, determining the cell identity and decoding the MIB and system information blocks (SIBs). This example demonstrates this process and decodes the MIB and SIB1, the first of the system information blocks. Decoding the MIB and SIB1 requires a comprehensive that is capable of demodulating and decoding the majority of the downlink channels and signals.

## Introduction

In order to communicate with the network the UE must obtain some basic system information. This is carried by the MIB and SIBs. The MIB carries the most essential system information:

- System bandwidth

- System frame number (SFN)

- Physical hybrid automatic repeat request (HARQ) indicator channel (PHICH) configuration

The MIB is carried on the broadcast channel (BCH) mapped into the physical broadcast channel (PBCH). This is transmitted with a fixed coding and modulation scheme and can be decoded after the initial cell search procedure. With the information obtained from the MIB the UE can now decode the control format indicator (CFI), which indicates the physical downlink control channel (PDCCH) length. This allows the PDCCH to be decoded, and searched for downlink control information (DCI) messages. A DCI message CRC masked with system information radio network temporary identifier (SI-RNTI) indicates that a SIB is carried in the same subframe. The SIBs are transmitted in the broadcast control channel (BCCH) logical channel. Generally, BCCH messages are carried on the downlink shared channel (DL-SCH) and transmitted on the physical downlink shared channel (PDSCH). The format and resource allocation of the PDSCH transmission is indicated by a DCI message on the PDCCH.

This example decodes the MIB and SystemInformationBlockType1 (SIB1). The latter is transmitted to specify the scheduling of other system information, along with aspects of the cell identity such as Public Land Mobile Network (PLMN) identity. Although SIB1 is transmitted in a fixed time schedule, the resource allocation of the PDSCH carrying SIB1 is dynamic and is indicated in an associated DCI message.

## Load and Process I/Q Waveform

In this example a single antenna I/Q capture of an eNodeB with two transmit antennas is used. The capture is performed at 15.36 Msamples/s which is sufficient to correctly sample all valid eNodeB bandwidths up to 10 MHz: 1.4 MHz, 3 MHz, 5 MHz, 10 MHz. The captured data is stored in the file eNodeBOutput.mat.

## Cell Search, Cyclic Prefix Length and Duplex Mode Detection

Call lteCellSearch to obtain the cell identity and timing offset offset to the first frame head. The cell search is repeated for each combination of cyclic prefix length and duplex mode, and the combination with the strongest correlation allows these parameters to be identified. A plot of the correlation between the received signal and the PSS/SSS for the detected cell identity is produced. The PSS is detected using time-domain correlation and the SSS is detected using frequency-domain correlation. Prior to SSS detection, frequency offset estimation/correction using cyclic prefix correlation is performed. The time-domain PSS detection is robust to small frequency offsets but larger offsets may degrade the PSS correlation.

<img width="250" alt="1" src="https://github.com/serttasbugrahan/lte-cell-search-mib-sib1-recovery/assets/140887398/c522f317-3f25-4c77-964d-da34cb4d5bd4">

## Frequency Offset Estimation and Correction

Prior to OFDM demodulation, any significant frequency offset must be removed. The frequency offset in the I/Q waveform is estimated and corrected using lteFrequencyOffset and lteFrequencyCorrect. The frequency offset is estimated by means of correlation of the cyclic prefix and therefore can estimate offsets up to +/- half the subcarrier spacing i.e. +/- 7.5kHz.

_Performing frequency offset estimation..._
_Frequency offset: -14.902Hz_

## OFDM Demodulation and Channel Estimation

The OFDM downsampled I/Q waveform is demodulated to produce a resource grid rxgrid. This is used to perform channel estimation. hest is the channel estimate, nest is an estimate of the noise (for MMSE equalization) and cec is the channel estimator configuration.

For channel estimation the example assumes 4 cell specific reference signals. This means that channel estimates to each receiver antenna from all possible cell-specific reference signal ports are available. The true number of cell-specific reference signal ports is not yet known. The channel estimation is only performed on the first subframe, i.e. using the first L OFDM symbols in rxgrid.

A conservative 13-by-9 pilot averaging window is used, in frequency and time, to reduce the impact of noise on pilot estimates during channel estimation.

## PBCH Demodulation, BCH Decoding, MIB Parsing

The MIB is now decoded along with the number of cell-specific reference signal ports transmitted as a mask on the BCH CRC. The function ltePBCHDecode establishes frame timing modulo 4 and returns this in the nfmod4 parameter. It also returns the MIB bits in vector mib and the true number of cell-specific reference signal ports which is assigned into enb.CellRefP at the output of this function call. If the number of cell-specific reference signal ports is decoded as enb.CellRefP=0, this indicates a failure to decode the BCH. The function lteMIB is used to parse the bit vector mib and add the relevant fields to the configuration structure enb. After MIB decoding, the detected bandwidth is present in enb.NDLRB.



