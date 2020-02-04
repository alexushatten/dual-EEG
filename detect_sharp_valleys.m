function ValE = detect_sharp_valleys(Signal)
% R Marquis @ FBM, November 2017
%
% very useful for pseudo-deterministic square-wave signals

Temp1 = diff(Signal);
    Temp2 = Temp1(1:end-1)>0; % UP (shifted)
    Temp1 = Temp1(2:end)<0; % DOWN (shifted)
    ValE = find((Temp1+Temp2)>1)+1; % UP AND, RIGHT AFTER, DOWN
%     figure; plot(STIMmkr(MkrChans(n),:));
%     hold on; plot(Temp3,STIMmkr(MkrChans(n),Temp3),'ro');

end
