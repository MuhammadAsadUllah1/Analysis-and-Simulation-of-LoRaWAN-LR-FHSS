% Paper title: Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario
% IEEE XPlore: https://ieeexplore.ieee.org/document/9653679
% Authors: Muhammad Asad Ullah, Konstantin Mikhaylov, Hirley Alves

% Cite this: M. A. Ullah, K. Mikhaylov and H. Alves, "Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario," in IEEE Wireless Communications Letters, doi: 10.1109/LWC.2021.3135984.

function [PS_DR8_analytical,H_N_Pro_succ,P_F] = DR8_analytical (N,pkct_p_h,Header_duration,F_duration,Last_fragment_duration,fragment_length,Header_N_DR8,Threshold,OBW_channels)

%% Probability of header success
% 1. Header collision with header
H_active_2TH = ((N*pkct_p_h*(Header_N_DR8))/3600)*2*Header_duration;
% 2. Header collision with payload data fragments
H_active_TF_TH = ((N*pkct_p_h*(fragment_length-Header_N_DR8-1-1))/3600)*(Header_duration+F_duration);
% 3. Header collision with last data fragment
H_active_TL_TH = ((N*pkct_p_h)/3600)*(Header_duration+Last_fragment_duration);
% Total colliding packet elements
H_active = H_active_2TH + H_active_TF_TH + H_active_TL_TH;


H_probability_succ=(1-(1/OBW_channels))^(H_active-1);
H_N_Pro_succ = 1-(1-H_probability_succ).^Header_N_DR8;

%% Fragments success probability, even if 1/3 (for DR8) of the fragments are lost

%% Payload data Fragment collisions
%1. Fragments collision with fragments
F_active_12_TF = ((N*pkct_p_h*(fragment_length-Header_N_DR8-1-1))/3600)*2*F_duration;
%2. Fragments collision with headers
F_active_12_TH = ((N*pkct_p_h*(Header_N_DR8))/3600)*(Header_duration+F_duration);
%3. Fragments collision with last data fragment
F_active_12_TL = ((N*pkct_p_h)/3600)*(F_duration+Last_fragment_duration);
% Total colliding packet elements
F_active_12 = F_active_12_TF + F_active_12_TH + F_active_12_TL;

Frag_probability_succ=(1-(1/OBW_channels))^(F_active_12-1);           % excluding last fragment
%% Last fragment collisions
% 1. Last fragment collision with last fragments
F_active_last_TL = ((N*pkct_p_h)/3600)*2*(Last_fragment_duration);
% 2. Last fragment collision with headers
F_active_last_TH = ((N*pkct_p_h*(Header_N_DR8))/3600)*(Last_fragment_duration+Header_duration);
% 3. Last fragment collision with payload data fragments
F_active_last_TF = ((N*pkct_p_h*(fragment_length-Header_N_DR8-1-1))/3600)*(F_duration+Last_fragment_duration);
% Total colliding packet elements
F_active_last = F_active_last_TH +F_active_last_TF+F_active_last_TL;

Frag_probability_succ_last=(1-(1/OBW_channels))^(F_active_last-1);    %Last fragment

Total_fragment = fragment_length-Header_N_DR8-1;
Frag_avg_suc = (sum(Frag_probability_succ*ones(1,(fragment_length-Header_N_DR8-2)))+Frag_probability_succ_last)/(Total_fragment);
%% Binomial probability
Fragment_success = binopdf(1:1:Total_fragment,Total_fragment,Frag_avg_suc);

%% Overall probability = Header success * fragment success
P_F = sum(Fragment_success(Threshold:1:Total_fragment));
PS_DR8_analytical = H_N_Pro_succ*sum(Fragment_success(Threshold:1:Total_fragment)); 
end