function [ts]=mylpfilt(ts,TR)
%______________________________
%function [ts]=mylpfilt(ts,TR);
%______________________________
% Pass time series and the TR in seconds.  Returns 'ts' low-pass filtered at 0.1Hz.
% if nargin<3
%     decision='n';
% end

[a,b]=lpfilt_creatorfcn(TR);
ts=filtfilt(b,a,ts);





%_________________________________________________________
function [a,b]=lpfilt_creatorfcn(TR)

fs=(1/TR);

Rp=3;%this is the allowable passband ripple in dB
Rs=60;%this is the required stopband attinuation in dB

Wp=0.1/(fs/2);
Ws=0.12/(fs/2);
[n,Wn]=cheb2ord(Wp,Ws,Rp,Rs);
[b,a]=cheby2(n,Rs,Wn);


% %this provides a figure displaying the lowpass filter using the current
% %parameters
% if decision=='d'
%     figure(1);
%     freqz(b,a,8192,fs);
%     title('Digital lowpass filter');
% end


