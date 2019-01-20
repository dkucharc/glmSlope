n      = 1000;
p      = 1000;

mFDR = 0;
sFDR = 0;
mPOW = 0;
sPOW = 0;

%n=p=1000 - 1000 powtórzeń około 30 minut

for q  = 0.1%[0.05,0.1,0.2]
    q_i    = 1-(1:p)*q/(2*p);
lambda = norminv(q_i,0,1);
lambda = lambda/4;
b_mean = zeros(p,1);
for p0 = [1 3 5 10 20 30 50]     %liczba niezerowych elementów wekt. b
FDR    = 0;
POWER  = 0;
b_est  = zeros(p,1);
%wekt = sign(randn(p,1));
for i=1:100

A      = randn(n,p);
A      = cat(2,orth(A-repmat(mean(A),size(A,1),1)),(1/sqrt(n)).*ones(n,1));
%A      = diag(ones(1,n));          %Alternatywnie macierz jednostkowa
p      = size(A,2);
mag    = 100;%25*sqrt(2*log(p));     %siła sygnału
b      = cat(1,repmat(mag,p0,1),zeros(p-p0,1));%.*wekt;

y      = binornd(1,1./(1+exp(-A*b)))*2-1;


[w,info] = logistic_slope(A,y,lambda);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 %Adlas(A/2,y,lambda);%


V     = sum(sign(w((p0+1):p).^2));
R     = max(cat(1,sum(sign(w.^2)),1));
FDR   = cat(1,FDR,V/R);

POWER = cat(1,POWER,sum(sign(w(1:p0).^2))/max(p0,1));

b_est = cat(2,b_est,w);
end
b_mean = cat(2,b_mean,mean(b_est(:,2:end),2));

FDR   = FDR(2:(i+1));
POWER = POWER(2:(i+1));

mFDR = cat(1,mFDR,mean(FDR));
sFDR = cat(1,sFDR,std(FDR)/sqrt(n));
mPOW = cat(1,mPOW,mean(POWER));
sPOW = cat(1,sPOW,std(POWER)/sqrt(n));
end
end
plot(b_mean)
mFDR(2:length(mFDR))
sFDR(2:length(sFDR))
mPOW(2:length(mPOW))
sPOW(2:length(sPOW))