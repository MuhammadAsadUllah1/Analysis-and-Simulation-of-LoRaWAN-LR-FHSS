% Paper title: Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario
% IEEE XPlore: https://ieeexplore.ieee.org/document/9653679
% Authors: Muhammad Asad Ullah, Konstantin Mikhaylov, Hirley Alves

% Cite this: M. A. Ullah, K. Mikhaylov and H. Alves, "Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario," in IEEE Wireless Communications Letters, doi: 10.1109/LWC.2021.3135984.

function [Payload_CRC_ToA_DR8,Payload_CRC_ToA_DR8_WH] = ToA_Packets_DR8(Payload,Header_ToA_DR8,M)    
   
    for PL=1:length(Payload)
    Payload_CRC_ToA_DR8(PL) = Header_ToA_DR8  + ceil((Payload(PL) + 2)/M)*(102/1000); 
    Payload_CRC_ToA_DR8_WH(PL) = ceil((Payload(PL) + 2)/M)*(102/1000); 
    end
end