function [Ttt,rd,si2,roi,d,n,pp,sit,si1,p,tt1,tt2] = SplitIntFunction_fin( T,Q ,xx)
%Calculates splitting intesity in frequency domain with error estimate from
%dominant frequencies
%   Uses the random windows chosen from Q and T

% Copyright 2024 F.Link and M.D.Long 
tau = xx;
Qn = Q';
Tn = T';
tn = tau.*(0:length(Qn)-1);
t1 = 0.1*round(length(Tn))*tau;
t2 = 0.9*round(length(Tn))*tau;
ind = 1:length(Qn);
nn=0;
sit2 = [];
pp2 = [];
for i = 1:500
    if i == 1
        tn1 = t1;
        tn2 = t2;
    else
        [tn1,tn2]  =  rand_wind(t1,t2,tn);
    end
    ioi = ind(tn>tn1&tn<tn2);
    Qt = Qn(ioi).*tukeywin(length(ioi),0.1);
    Tt = Tn(ioi).*tukeywin(length(ioi),0.1);

    % Equation from Chevrot, 2000
    rd = zeros(length(Tt)-1,1);
    for put = 1:length(Tt)-1
        rd(put)= (Qt(put+1)-Qt(put))/xx;
    end
    Ttt = Tt;
    Ttt(end) = [];
    Pop= dot(Ttt,rd);
    sn=sign(Pop);
    NPop =norm(rd)^2;
    si2(i) = -2.*Pop./NPop;

    % Frequency dependent SI

    L = length(ioi);
    L2 = 2.^(nextpow2(L)+6);

    Qt2 = zeros(L2,1);
    Qt2(1:L) = Qt;
    Tt2 = zeros(L2,1);
    Tt2(1:L) = Tt;
    om = (0:fix(L2/2)-1)./L2./tau.*2.*pi;

    rr = fft(Qt2);
    d = rr(1:fix(L2/2)).*sqrt(-1).*om';
        roi = abs(rr(1:fix(L2/2)));
    o1 = 2*pi/6;
    o2 = 2*pi/8;
    pn = sum(om>o2&om<o1);

    % get dominant frequency

    TT = fft(Tt2);
    n = TT(1:fix(L2/2));

    huh0 = abs(conj(d).*d);
    hm0 = real(conj(n).*d);

    hm = movmean(hm0,[pn/2 pn/2]);
    huh = movmean(huh0,[pn/2 pn/2]);


    sit = -2.*hm./huh;
    pp = 2*pi./om;

    [~,I] = max(abs(roi));
    si1(i) = sit(I);

    p(i) = pp(I);

    tt1(i) = tn1;
    tt2(i) = tn2;

end


end