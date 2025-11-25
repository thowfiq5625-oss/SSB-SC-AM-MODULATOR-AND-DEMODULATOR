# SSB-SC-AM-MODULATOR-AND-DEMODULATOR
## AIM

To write a program to perform SSBSC modulation and demodulation using SCI LAB and study its spectral characteristics
## EQUIPMENTS REQUIRED

•	Computer with i3 Processor

•	SCI LAB

Note: Keep all the switch faults in off position

## ALGORITHM
  1.	Define Parameters:
     
    •	Fs: Sampling frequency.
    •	T: Duration of the signal.
    •	Fc: Carrier frequency.
    •	Fm: Frequency of the message signal.
    •	Amplitude: Maximum amplitude of the message signal.
    
  2.	Generate Signals:
     
    •	Message Signal: The baseband signal that will be modulated.
    •	Carrier Signal: A high-frequency signal used for modulation.
    •	Analytic Signal: Constructed using the Hilbert transform to get the in-phase and quadrature components.

  3.	SSBSC Modulation:
     
    •	Modulated Signal: Create the SSBSC signal using the in-phase and quadrature components, modulated by the carrier.

  4.	SSBSC Demodulation:
     
    •	Mixing: Multiply the SSBSC signal with the carrier to retrieve the message signal.
    •	Low-pass Filtering: Apply a low-pass filter to remove high-frequency components and recover the original message signal.

  5.	Visualization:
     
    Plot the message signal, carrier signal, SSBSC modulated signal, and the recovered signal after demodulation.

## PROCEDURE

  •	Refer Algorithms and write code for the experiment.
  
  •	Open SCILAB in System
  
  •	Type your code in New Editor
  
  •	Save the file
   
  •	Execute the code
  
  •	If any Error, correct it in code and execute again
  
  •	Verify the generated waveform using Tabulation and Model Waveform
## MODEL GRAPH
<img width="747" height="338" alt="image" src="https://github.com/user-attachments/assets/dd117c5d-ee32-47c7-946c-df6180b0d33f" />

## PROGRAM
clc;
clear;
close;

Am = 10.5;

Ac = 21; 

fm = 1033; 

fc = 10330; 

fs = 1033000; 

t = 0:1/fs:2/fm;

m1 = Am * cos(2 * %pi * fm * t);

c1 = Ac * cos(2 * %pi * fc * t);

s1 = c1 .* m1;

m2 = Am * cos(1.57 - (2 * %pi * fm * t));

c2 = Ac * cos(1.57 - (2 * %pi * fc * t));

subplot(5,1,1);

plot(t, m1);

xtitle("Message Signal", "Time (s)", "Amplitude");

subplot(5,1,2);

plot(t, c1);

xtitle("Carrier Signal", "Time (s)", "Amplitude");

s2 = c2 .* m2;

lsb = s1 + s2;

usb = s1 - s2;

subplot(5,1,3);

plot(t, lsb);

xtitle("Lower Sideband (SSB-LSB)", "Time (s)", "Amplitude");

subplot(5,1,4);

plot(t, usb);

xtitle("Upper Sideband (SSB-USB)", "Time (s)", "Amplitude");

ssb_signal = usb;

demod = ssb_signal .* cos(2 * %pi * fc * t);

demod = demod / (Ac / 2); 

fc_cut = 2 * fm;

N = 101;

M = (N-1)/2;

n = 0:N-1;

function y = mysinc(x)

    y = zeros(x);
    
    for i = 1:length(x)
    
        if x(i) == 0 then
          
            y(i) = 1;
        
        else
            y(i) = sin(%pi * x(i)) / (%pi * x(i));
        end
    end
endfunction

h = (2 * fc_cut / fs) * mysinc(2 * fc_cut * (n - M) / fs);

w = 0.54 - 0.46 * cos(2 * %pi * n / (N - 1));

h = h .* w;

h = h / sum(h);

recovered = filter(h, 1, demod);

subplot(5,1,5);

plot(t, recovered, 'r');

xtitle("Recovered Message Signal (After LPF)", "Time (s)", "Amplitude");

## TABULATION
![WhatsApp Image 2025-11-25 at 16 11 31_5bc46e78](https://github.com/user-attachments/assets/c9b8cd90-de3f-40b7-998c-e6e36c95ebc4)

## OUTPUT
![WhatsApp Image 2025-11-19 at 17 45 51_5fb9ca6b](https://github.com/user-attachments/assets/089ba94b-39a3-4e04-8ac7-12566f594212)

## RESULT
Thus, the SSB_SC_AM Modulation and Demodulation is experimentally done and the output is verified.


