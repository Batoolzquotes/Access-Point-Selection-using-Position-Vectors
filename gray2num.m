function a=gray2num(g)
for i=1:length(g)
    t=0;
    for j=1:i
        t=bin2dec(g(j))+t;
    end
    d(i)=dec2bin(mod(t,2));
end
a=bin2dec(d);
return 