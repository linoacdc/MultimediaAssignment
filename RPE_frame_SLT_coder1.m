
t = 1.25*10^-4:1.25*10^-4:0.02;
s0=sin(100000*t)+sin(5000*t);
PrevFrmSTResd =zeros(160,1);
[LARc,Nc,bc,CurrFrmExFull,CurrFrmSTResd] = RPE_frame_SLT_coder(s0,PrevFrmSTResd);
LARc,Nc,bc,CurrFrmExFull,CurrFrmSTResd;
[s,CurrFrmSTResd]= RPE_frame_SLT_decoder(LARc,Nc,bc,CurrFrmExFull,PrevFrmSTResd);
s,CurrFrmSTResd

function [LARc,Nc,bc,CurrFrmExFull,CurrFrmSTResd] = RPE_frame_SLT_coder(s0,PrevFrmSTResd)


%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%
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




%%Calculate s' and d
%s^(n)=Sum(ak*s(n-k))


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

%Calculate difference between s and sDec, I dont use this
 d1=s-sDec;
 
%Calculate difference based on formula of the standard given, this is 
%what I use
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
 
%Short term analysis residual
d2=di(:,9);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LONG TERM ANALYSIS

%Initialize values
Prevd = zeros(120,1);
d = zeros(40,1);
N = zeros(4,1);
b = zeros(4,1);
dRecon = zeros(160,1);
dReconCurrent = zeros(40,1);
dEst = zeros(40,1);
e = zeros(160,1);
QLB = [0.1, 0.35, 0.65, 1.0];
bc=zeros(4,1);
Nc=zeros(4,1);
bdec = zeros(4,1);
Ndec = zeros(4,1);

%Calculate N,b,d',e for each subframe
for j = 0:3
    
    count = 0;
    d = d2(1+j*40:j*40+40);

    %Setting up Prevd and d for each subframe to use for RPE_subframe
    %function
    for i = 3:-1:1
        if (j*40 + 1 - i*40 <= 0)
            PrevFrmSTResd(j*40 + 1 - i*40 +160: j*40 - i*40 + 40 + 160)=dReconCurrent(1:40);
            Prevd(count*40 + 1:count*40 + 40) = dReconCurrent(1:40);
        elseif(j*40 + 1 - i*40 > 0)
            Prevd(count*40 + 1:count*40 + 40) = dRecon(j*40 + 1 - i*40:j*40 - i*40 +40);
        end
        count=count + 1;
    end
    
    %Use function to find N, b
    [N(j+1),b(j+1)] = RPE_subframe_LTE(d, Prevd );
    
    %Encode N, b
    Nc(j+1) = N(j+1);
        if(b(j+1) <= 0.2)
            bc(j+1) = 0;
        elseif(0.2 < b(j+1) <= 0.5)
            bc(j+1) = 1;
        elseif(0.5 < b(j+1) <= 0.8)
            bc(j+1) = 2;
        elseif(0.8 < b(j))
            bc(j+1) = 3;
        end
    
    %Decode N, b
    Ndec(j+1) = N(j+1);
    bdec(j+1) = QLB(bc(j+1)+1);
    
    %Calculate d'' like in the assignment presentation d'' = bc*d'(n-Nc)
    dEst(1:40) = bc(j+1) * Prevd(120 - Nc(j+1) + 1: 120 - Nc(j+1) +40);
   
    %Calculate e(n) = d(n) - d''(n)
    e(j*40 +1: j*40 + 40) = d2(j*40 +1: j*40 + 40) - dEst(1:40);
    
    %Calculate reconstructed d' = e(n)' + bdec * d'(n-Ndec), e(n) = e'(n)
    %for this assignment
    dRecon(j*40 +1: j*40 + 40) = e(j*40 +1: j*40 + 40) + bdec(j+1) * Prevd(120 - Ndec(j+1) + 1: 120 - Ndec(j+1) +40);
    
end
 
CurrFrmSTResd = dRecon;
CurrFrmExFull = e;

end

function [s0, CurrFrmSTResd] = RPE_frame_SLT_decoder(LARc,Nc,bc,CurrFrmExFull, PrevFrmSTResd)

%%%Initialize some values
bdec = zeros(4,1);
Ndec = zeros(4,1);
e = (CurrFrmExFull);
QLB = [0.1, 0.35, 0.65, 1.0];
d = zeros(160,1);
Prevd = PrevFrmSTResd(41:160);

%Calculate d' per subframe from d?(n) = e?(n) + b?d?(n-N')
for j = 0:3
   
    %Decode N,b
    Ndec(j+1) = Nc(j+1);
    bdec(j+1) = QLB(bc(j+1)+1);
    
    %Calculate d for subframe
    d(40*j + 1 :40*j + 40) = e(40*j + 1 :40*j + 40) + bdec(j+1) * Prevd(120 - Ndec(j+1) + 1: 120 - Ndec(j+1) +40);
    
    %Change Prevd to match the next subframe
    Prevd(1:80) = Prevd(41:120);
    Prevd(81:120) = d(40*j + 1 :40*j + 40);
end

CurrFrmSTResd = d;

%Calculate reflection coefficients from LARc
A =  [20, 20, 20, 20, 13.637, 15, 8.334, 8.824];
B =  [0, 0, 4, -5, 0.184, -3.5, -0.666, -2.235];
MinLar = [-32, -32, -16, -16, -8, -8, -4, -4];
MaxLar = [31, 31, 15, 15, 7 ,7, 3, 3];


%Decode LARc into LARdec
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

%Calculate a(not needed here)
aDec = rc2poly(rDec);

%Initialize values for short term filter
s = zeros(160,1);
sr = zeros(160,9);
v = zeros(160,9);


%Use d through the filter to get s
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


s0 = zeros(160,1);
s0(1)=s(1);

%Postprocessing
for k = 2:160
    s0(k) = s(k) + 28180*2^(-15) * s(k-1);
end


end

function [N,b] = RPE_subframe_LTE(d, Prevd )

sumLambda = zeros(81,1);

%Calculating all cross-correlations
for lambda = 40:120
    for i = 1:40
        sumLambda(lambda-39) = sumLambda(lambda-39) + Prevd(120-lambda+i)*d(i);
end
end

max= sumLambda(1);
N=40;

%Find the largest and store the index for that ( add 39 so it matches with
%N)
for i =2:81
    if (sumLambda(i)>=max)
        max = sumLambda(i);
        N=i+39;
    end
end

sum1=0;
sum2=0;

%Calculate the two sums needed to calculate b
for i = 1:40
    sum1 = sum1 + d(i)*Prevd(120 + i - N);
    sum2 = sum2 + (Prevd(120 + i - N)^2);
end
    
b = sum1/sum2;


end

function zrounded = Nint(z)
    
    %zrounded = int16(z+sign(z)*0.5);
    zrounded = round(z);

end
