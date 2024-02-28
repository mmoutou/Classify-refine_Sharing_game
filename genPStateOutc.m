% just to check weighed sum of outcome probabilities :)

AH = 100; BH = 100;
Ns = 10000;
Ao=[30 2]; Bo=[10 4];
No = 100;
gpo = zeros(Ns,1);
for k = 1:Ns
    pH = betarnd(AH,BH,1);
    % pS = betarnd(AS,BS,1);
    s = pBinSample([1-pH, pH],1); 
    gpo(k) = betarnd(Ao(s),Bo(s),1);
end
hist(gpo,30)