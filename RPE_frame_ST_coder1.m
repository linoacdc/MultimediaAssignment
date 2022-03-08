
t = 1.25*10^-4:1.25*10^-4:0.02;
s0=sin(1000*t)+0.5*sin(5000*t)+0.7*sin(4000*t)+sin(200*t);
figure(1);
plot(s0);
title('s0');

[LARci,CurrFrmSTResd] = RPE_frame_ST_coder(s0);
LARci,CurrFrmSTResd
[s]= RPE_frame_ST_decoder(LARci,CurrFrmSTResd);
figure(2);
plot(s);
title('s');

function [LARc,CurrFrmSTResd] = RPE_frame_ST_coder(s0)

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREPROCESSING%
alpha = 32735*2^-15;
beta = 28180*2^-15;
s=zeros(1,160);
s(1)=s0(1);

for k = 2:160
    
    s(k) = s0(k) - s0(k-1) + alpha*s(k-1);
    
end
for k = 2:160
    
    s(k) = s(k) - beta*s(k-1);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SHORT TERM ANALYSIS

%Calculate estimation of autocorrelation values r=ACF
 ACF=zeros(1,9);
 w=zeros(1,8);
 
 for k = 1:9
    for i = k+1:160
        ACF(k)=ACF(k) + s0(i)*s0(i-k);
    end
 end

%Calculating the specific autocorrelation matrix R and r
for i = 1:8
    for j = 1:8
        R(i,j)=ACF(abs(i-j)+1);
        
    end
    r(i)=ACF(i+1);
end


%Solve linear system to find the optimal coefficients
w=linsolve(R,transpose(r));


% Transform coefficients to reflection coefficients w->kr

%Create vector for poly2rc
a=zeros(9,1);
a(1)=1;
a(2:9)=-w(1:8);


%Use poly2rc to get reflection coefficients kr
kr=poly2rc(a);

%Check and make sure max values are within (-1,1)
gt1 = kr > 1;
ltm1 = kr < -1;
kr(gt1)=1;
kr(ltm1)=-1;

%Create LAR from reflection coefficients
for i=1:8
    if abs(kr(i)) >= 0.950 && abs(kr(i)) <= 1
        LAR(i) = sign(kr(i)) * (8*abs(kr(i))-6.375);
    elseif abs(kr(i)) >= 0.675 && abs(kr(i))< 0.950
        LAR(i) = sign(kr(i)) * (2*abs(kr(i))-0.675);
    elseif abs(kr(i)) < 0.675
        LAR(i)=kr(i);
    end
end
% LAR

%Could use this but standard uses the approximation above
% for i=1:8
%     LAR(i)=log10((1+kr(i))/(1-kr(i)));
% end


%Quantize and integet code LAR into LARc
A =  [20, 20, 20, 20, 13.637, 15, 8.334, 8.824];
B =  [0, 0, 4, -5, 0.184, -3.5, -0.666, -2.235];
MinLar = [-32, -32, -16, -16, -8, -8, -4, -4];
MaxLar = [31, 31, 15, 15, 7 ,7, 3, 3];

for i = 1:8
    LARc(i) = Nint(A(i)*LAR(i) + B(i));
    if LARc(i) > MaxLar(i)
        LARc(i)= MaxLar(i);
    elseif LARc(i) < MinLar(i)
        LARc(i) = MinLar(i);
    end
end

%Calculate decoded LAR (LARc->LARdec)
for i=1:8
    LARdec(i)=( (LARc(i)) - B(i) ) / A(i);
end


%Calculate decoded reflection coefficients
for i=1:8
    if abs(LARdec(i)) <= 1.625 && abs(LARdec(i)) >= 1.225
        rDec(i)=sign(LARdec(i)) * (0.125*abs(LARdec(i)) + 0.796875);
    elseif abs(LARdec(i)) < 1.225 && abs(LARdec(i)) >= 0.675
        rDec(i)=sign(LARdec(i)) * (0.5*abs(LARdec(i)) + 0.337500);
    elseif abs(LARdec(i)) < 0.675
        rDec(i)=LARdec(i);
    end
end

%Calculate coefficients from reflection coefficients
aDec = rc2poly(rDec);
sDec = zeros(1,160);

%Calcualte decoded signal from formula
for n=1:160
    for k=1:8
        if n>=k
            sDec(n)=sDec(n)+ aDec(k) * s(n-k+1);
        end
    end
end

%Calculate difference between s and sDec

 d1=s-sDec;

ui=zeros(160,9);
di = zeros(160,9);

for i = 1:9
    di(1,i)=s(1);
    ui(1,i)=s(1);        
end

 for k=2:160
     di(k,1)=s(k);
     ui(k,1)=s(k);
     
     for i=2:9
             di(k,i) = di(k,i-1) + rDec(i-1)*ui(k-1,i-1);
             ui(k,i) = ui(k-1,i-1) + rDec(i-1) * di(k,i-1);
     end
 end
d2=di(:,9);
CurrFrmSTResd=d2;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s0] = RPE_frame_ST_decoder(LARc,CurrFrmSTResd)

%Decode LAR
A =  [20, 20, 20, 20, 13.637, 15, 8.334, 8.824];
B =  [0, 0, 4, -5, 0.184, -3.5, -0.666, -2.235];
MinLar = [-32, -32, -16, -16, -8, -8, -4, -4];
MaxLar = [31, 31, 15, 15, 7 ,7, 3, 3];

for i=1:8
    LARdec(i)=( (LARc(i)) - B(i) ) / A(i);
end

%Calculate rDec from LARdec
for i=1:8
    if abs(LARdec(i)) <= 1.625 && abs(LARdec(i)) >= 1.225
        rDec(i)=sign(LARdec(i)) * (0.125*abs(LARdec(i)) + 0.796875);
    elseif abs(LARdec(i)) < 1.225 && abs(LARdec(i)) >= 0.675
        rDec(i)=sign(LARdec(i)) * (0.5*abs(LARdec(i)) + 0.337500);
    elseif abs(LARdec(i)) < 0.675
        rDec(i)=LARdec(i);
    end
end

%Calculate aDec
aDec = rc2poly(rDec)

s = zeros(160,1);
sr = zeros(160,9);
v = zeros(160,9);

d = transpose(CurrFrmSTResd);
%Calculate s from rDec and current d
for i = 1:9
    sr(1,i)=d(1);  
end
s(1)=sr(1,9);

for k = 2:160

    sr(k,1) = d(k);    
    
    for i = 2:9
        
        sr(k,i) = sr(k,i-1) - rDec(10-i)*v(k-1,10-i);
        v(k,11-i) = v(k,10-i) + rDec(10-i)*sr(k,i);
    end
    s(k) = sr(k,9);
    v(1,k) = sr(k,9);
end

%Postprocessing
s0 = zeros(160,1);

s0(1)=s(1);
for k = 2:160
    s0(k) = s(k) + 28180*2^(-15) * s(k-1);

end
end

function zrounded = Nint(z)
    zrounded = round(z);
end
