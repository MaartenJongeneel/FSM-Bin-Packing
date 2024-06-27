function [H_filt,V,dV] = FilteringSE3(H, t, fwp, fwo, p)
%UNTITLED Summary of this function goes here
%   Detailed explanation go here
% https://nl.mathworks.com/help/signal/ug/take-derivatives-of-a-signal.html
       
    fc = round(1/mean(diff(t)));

    [R,w,dw,~] = sgolayfiltSO3(H(1:3,1:3,:),p, fwo, fc);
    
    R = cat(3,repmat(R(:,:,1),1,1,fwo+1),R);
    R = cat(3,R,repmat(R(:,:,end),1,1,fwo));
    w = [zeros(3,fwo+1) w zeros(3,fwo)];
    dw = [zeros(3,fwo+1) dw zeros(3,fwo)];
    
    o = smoothdata(permute(H(1:3,4,:), [3,1,2]),'sgolay',fwp);

    Nf = 50; 
    Fpass = 20; 
    Fstop = 60;
    d = designfilt('differentiatorfir','FilterOrder',Nf, ...
        'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
        'SampleRate',fc);
    
    dt = diff(t);
    dt(end+1) = mean(diff(t));
    do = filter(d,o)./dt;
    
    % The filtered signal is delayed. 
    % Use grpdelay to determine that the delay is half the filter order. 
    % Compensate for it by discarding samples.
    % The output also includes a transient whose length equals the filter 
    % order, or twice the group delay. delay samples were discarded above. 
    % Discard delay more to eliminate the transient. (2*delay)
    delay = mean(grpdelay(d));
    do(1:2*delay,:) = zeros(2*delay,3);
    ddo = filter(d,do)./dt;
    ddo(1:4*delay,:) = zeros(4*delay,3);

    V = [do'; w];
    dV = [ddo'; dw];
    
    H_filt = zeros(4,4,length(H));
    for idx=1:length(H)
        H_filt(:,:,idx) = [R(:,:,idx), o(idx,:)'; 0 0 0 1];       
    end

end

