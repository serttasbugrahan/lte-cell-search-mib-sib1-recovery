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

![GitHub Logo](C:\Users\BUGRAHAN\Desktop\1.png)
