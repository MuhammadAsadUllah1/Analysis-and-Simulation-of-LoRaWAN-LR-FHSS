# Analysis-and-Simulation-of-LoRaWAN-LR-FHSS
# Analysis and Simulation of LoRaWAN LR-FHSSfor Direct-to-Satellite Scenario

This is a code is related to the following article:

M. Asad Ullah, K. Mikhaylov and H. Alves, "[Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario](https://ieeexplore.ieee.org/document/9653679)," in IEEE Wireless Communications Letters, doi: 10.1109/LWC.2021.3135984.

In this code package, ***LR_FHSS_DR8_Capture_main_file.m*** is main script file which defines the simulation parameters. The function ***Satellite_Geometry(H,E)*** calculates the distance from node to satellite as function of elevation angle. Similarly, function ***ToA_Packets_DR8*** calculate the time-on-air (ToA) for Data Rate DR8. The rest of code follows the LR-FHSS functionalities to investigate the Packet Success Probability. Note that the current script files are for DR8 but it can be easily modify to other DRs 9-11 following the comments given in the code.

# Abstract
Long Range-Frequency Hopping Spread Spectrum (LR-FHSS) has been recently introduced into the LoRaWAN protocol specification to increase network capacity and collision robustness, and enable direct connectivity between machine devices and the Low Earth Orbit (LEO) satellites. In this letter, we first construct the analytical and simulation models for packet delivery over LR-FHSS from ground nodes to a LEO satellite, and then use the developed analytic and simulation models to generate the numerical results. Our results reveal the potential feasibility of large-scale networks, demonstrate some trade-offs between the two new LR-FHSS-based data rates for the EU region, and reveal the key reasons for packet losses.

# What is LoRaWAN LR-FHSS modulation?

With the motivation to overcome the connectivity gaps in remote areas, LoRa Alliance has introduced LR-FHSS data rates into the LoRaWAN protocol in the last quarter of 2020. The new data rates exploit the frequency hopping and offers high robustness against the co-channel interference through increasing the number of physical channels with redundant physical headers. In this context, the key point of our work is to provide the basis allowing to analyse the scalability of LR-FHSS. In the follow-up studies, we aim to approach the the associated implementation challenges.

The novel LR-FHSS scheme exploits frequency hopping and offers high robustness against the interfering through increasing the number of physical channels with redundant physical headers. Specifically, LR-FHSS uses the following mechanisms:

- Redundancy and coding: Unlike LoRa, the LR-FHSS device transmits severalreplicas of headers (**N = 1, . . . , 4**), where the number of replication is defined
by the date rate setting. To successfully decode a packet, a gateway should receive at least one of the N transmitted headers. For the payload data fragments,
DR8/DR10 and DR9/DR11 imply coding with the rate equivalent to **CR = 1/3** and **CR = 2/3** which improves the the gateway’s ability to correctly demodulate
the radio packets due to lower coding rate.

- Frequency hopping: For the uplink communication in the EU region, an LR-FHSS based end-device sends each packet fragment (i.e., the header or a data
fragment) on another randomly-picked frequency channel. In a single Operating Channel Width (OCW) channel, LR-FHSS specifications define **280** and **688** Occupied Band Width (OBW) physical carries (with a bandwidth of 488 Hz) for
DR8/DR9 and DR10/DR11, respectively. On contrary, LoRa modulation carry a complete transmission in a single channel which causes higher co-channel interference probability.


Due to all these promising features, LR-FHSS is a prominent technology enabling mMTC-over-satellite network to fill the connectivity gaps in hard-to-reach areas, which
lacks the terrestrial networks infrastructure e.g., offshore wind farms, and vessels monitoring in Arctic.

# LoRaWAN LR-FHSS Simulator

The LR-FHSS simulator operation can be subdivided into the three main phases:

- ***Time-frequency scheduling:*** For the uplink communication in the EU region, an
LR-FHSS based end-device initiates the transmission by sending **N** copies of the
0.233 ms long header on randomly selected carrier frequencies. Next, it breaks
the **L bytes** physical payload into several fragments of 50 ms and implies random
selection of the sub-channel for each transmission following pseudo-random generator (PRNG) for frequency hopping pattern. Note that the DR8/DR10 and
DR9/DR11 imply coding with the rate equivalent to **CR = 1/3** and **CR = 2/3**,
respectively, which improves the gateway’s ability to correctly demodulate the
hundreds of radio packets. Our simulator models this procedure in accordance
with the LoRaWAN specification documents.´
- ***Packet elements collision analysis:*** A node of interest transmits a single packet
element over one of the **C = 280** Occupied Band Width (OBW) channels at time
instant **τi,j** . At the same time, there might be other devices trying to utilize the
same OBW channel, which the target node is using for transmitting the target
header. Considering random channel selection and ALOHA-like media access,
probability of accessing the same OBW frequency channel is **1/C**. If a collision occurs (i.e., there is an overlap in time and frequency channel for two or more
elements), our simulator calculates and examines the received power of the target
packet element as well as the interfering packets elements at the gateway. We
imply that the gateway can successfully recover the target packet element from
a collision benefiting from the capture effect if the desired packet element is at
least δ (considered as **6 dB** in our work) dB stronger than the interfering packets.
By the end of this phase, for each packet element (i.e., header and data fragment)
the simulator has information of whether this packet element has been received
by the gateway or not.

- ***Packet delivery Ratio (PDR):*** To successfully decode a packet, it is essential that
the gateway receives at least one of the N transmitted headers and the number of received data fragments should be higher than the **pre-defined reception
threshold (γ)**, which depends on the coding rate. At the last stage of operation,
the simulator analyzes the number of headers and fragments received and makes
the final decision of whether the whole packet can be decoded. The results are
further logged

![image](https://user-images.githubusercontent.com/67004968/155481649-468d7b08-daf1-44c8-bb3b-ca08c5bf1e48.png)

# Acknowledgements
This work was financially supported byA cademy of Finland 6Genesis Flagship (Grant Number: 318927) and Academy of Finland MRATSafeDrone (Grant Number: 341111).

# Further relevant works
- M. Asad Ullah, K. Mikhaylov and H. Alves, "[Massive Machine-Type Communication and Satellite Integration for Remote Areas](https://ieeexplore.ieee.org/document/9535456)," in IEEE Wireless Communications, vol. 28, no. 4, pp. 74-80, August 2021, doi: 10.1109/MWC.100.2000477.
- M. Asad Ullah, K. Mikhaylov and H. Alves, "[Enabling mMTC in Remote Areas: LoRaWAN and LEO Satellite Integration for Offshore Wind Farm Monitoring](https://ieeexplore.ieee.org/document/9537682)," in IEEE Transactions on Industrial Informatics, vol. 18, no. 6, pp. 3744-3753, June 2022, doi: 10.1109/TII.2021.3112386.
# To Cite
M. Asad Ullah, K. Mikhaylov and H. Alves, "[Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario](https://ieeexplore.ieee.org/document/9653679)," in IEEE Wireless Communications Letters, doi: 10.1109/LWC.2021.3135984.
