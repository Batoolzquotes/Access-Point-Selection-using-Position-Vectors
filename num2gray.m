function c = num2gray(a,L)
b=dec2bin(a,L);
 c(1)=b(1);
for i=1:L-1
    m=bin2dec(b(i));
    n=bin2dec(b(i+1));
    c(i+1)=dec2bin(xor(m,n));
end

