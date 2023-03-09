function [SST2,o,oo,tt,phi22p] = sstn_test3(s,gamma,sigma,ft,bt)
n = length(s);
nv = log2(n);
if mod(nv,1)~=0
    warning('The signal is not a power of two, truncation to the next power');
    s = s(1:2^floor(nv));
end
n = length(s);
s = s(:);
if nargin<5
   ft = 1:n/2;
   bt = 1:n;
end
nb = length(bt);
neta = length(ft);
sz=zeros(n,1);
sleft = flipud(conj(sz(2:n/2+1)));
sright = flipud(sz(end-n/2:end-1));
x = [sleft; s ; sright];
clear xleft xright;
t = -0.5:1/n:0.5-1/n;t=t';
g =  1/sigma*exp(-pi/sigma^2*t.^2);
gp = -2*pi/sigma^2*t .* g; % g'
gpp = (-2*pi/sigma^2+4*pi^2/sigma^4*t.^2) .* g; % g''
SST2 = zeros(neta,nb);
o = zeros(neta,nb);
tt = zeros(neta,nb);
oo = zeros(neta,nb);
vg = zeros(neta,7);
vgp = zeros(neta,5);
Y = zeros(neta,4,4);
for b=1:nb	
    for i = 0:7
        tmp = (fft(x(bt(b):bt(b)+n-1).*(t.^i).*g))/n;
        vg(:,i+1) = tmp(ft);
    end       
    for i = 0:5
        tmp = fft(x(bt(b):bt(b)+n-1).*(t.^i).*gp)/n;
        vgp(:,i+1) = tmp(ft);
    end  
    tt(:,b) = vg(:,2)./vg(:,1);
    for i = 1:7
        for j = 1:7
            if i>=j
                Y(:,i,j) = vg(:,1).*vg(:,i+1) - vg(:,j).*vg(:,i-j+2);
            end
        end
    end      
    W2 = 1/2/1i/pi*(vg(:,1).^2+vg(:,1).*vgp(:,2)-vg(:,2).*vgp(:,1));
    o(:,b) = (ft-1)'-real(vgp(:,1)/2/1i/pi./vg(:,1));
    p(:,b) = W2./Y(:,2,2);
    oo(:,b) = o(:,b) + real(p(:,b).*tt(:,b));
    nn = exp(1i*pi*(ft-1)');
    aa(:,b) = vg(:,1).* nn;
    for eta=1:neta
        if abs(aa(eta,b))>gamma  
            k = 1+round(oo(eta,b));
             if k>=1 && k<=neta
               SST2(k,b)  = SST2(k,b) + aa(eta,b);
             end
        end
    end
end