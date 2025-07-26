M = % window length 
N = % length of the time series
t = (1:N); 
rep = %number of reconstructed components used for forecast  
fop =%number of points to be forecasted 
%Reading dew point temperature file 
dop = xlsread('inputfile_path');      
%Reading relative humidity file
vi = xlsread('inputfile_path2');
%finding mean of dew point temperature and relative humidity
meandop = mean(dop);
meanrh = mean(vi);
%finding standard deviation of dew point temperature and relative humidity
stddop = std(dop);
stdrh = std(vi);

%normalizing each of the time series by subrtracting mean from the times series and dividing by the times series by standard deviation 
dop = dop-mean(dop); % remove mean value
vi = vi-mean(vi);
dop = dop/std(dop);  % normalize to std=1
vi = vi/std(vi);
%combining both the time series into a single array 
X = [dop vi]; % multivariate time series

%plotting both the time series individually 
figure(1);
clf;
set(1,'name','Time series of dp and rh');
subplot(1,2,1);
plot(dop, 'r-');
title('Time series dp');
subplot(1,2,2);
plot(vi, 'r-');
title('Time series rh');

%calculating the covariance of both the time series by creating two zero matrices 
Y1=zeros(N-M+1,M);
Y2=zeros(N-M+1,M);
for m=1:M                 % create time-delayed embedding of X
  Y1(:,m) = dop((1:N-M+1)+m-1);
  Y2(:,m) = vi((1:N-M+1)+m-1);
end
Y = [Y1 Y2];
Cemb=Y'*Y / (N-M+1);

%plotting the covariance matrix in a color bar graph 
figure(2);
imagesc(Cemb);
axis square
set(gca,'clim',[-1 1]);
colorbar

C=Cemb;
xlswrite('outputfile',C,'covariance')

%calculating the eigen value and eigen vector of the covariance matrix 
[RHO,LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);      % extract the diagonal
[LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues and eigenvectors
RHO = RHO(:,ind);             
%calculating cumulative variance 
LAMBDA_RESULT=LAMBDA;
TOTAL_SUM_OF_EIGEN_VALUES=sum(LAMBDA_RESULT,'all');
cum_var=LAMBDA_RESULT/TOTAL_SUM_OF_EIGEN_VALUES;
cum_var_percentage=cum_var*100;
%disp(cum_var_percentage)
xlswrite('outputfile',cum_var_percentage,'cummulative covariance')
xlswrite('outputfile',LAMBDA,'eigen value')
xlswrite('outputfile',RHO,'eigen vector')

%plotting eigen value and eigen vector 
figure(3);
clf;
set(gcf,'name','Eigenvectors RHO and eigenvalues LAMBDA')
subplot(3,1,1);
plot(LAMBDA,'o-');
subplot(3,1,2);
plot(RHO(:,1:2), '-');
legend('1', '2');
subplot(3,1,3);
plot(RHO(:,3:4), '-');
legend('3', '4');

%calculating the principal components 
PC = Y*RHO;
%plotting the principal components 
figure(4);
set(gcf,'name','Principal components PCs')
clf;
for m=1:4
  subplot(4,1,m);
  plot(PC(:,m),'k-');
  ylabel(sprintf('PC %d',m));
  ylim([-10 10]);
end
xlswrite('outputfile',PC,'PC')

%calculating the reconstructed components RC1 and RC2 
RC1=zeros(N,2*M);
RC2=zeros(N,2*M);
for m=1:2*M
  buf1=PC(:,m)*RHO(1:M,m)'; % invert projection - first channel
  buf1=buf1(end:-1:1,:);

  buf2=PC(:,m)*RHO(M+1:end,m)'; % invert projection - second channel
  buf2=buf2(end:-1:1,:);

  for n=1:N % anti-diagonal averaging
    RC1(n,m)=mean( diag(buf1,-(N-M+1)+n) );
    RC2(n,m)=mean( diag(buf2,-(N-M+1)+n) );
  end
end

xlswrite('outputfile',RC1,'RC1')
xlswrite('outputfile',RC2,'RC2')
figure(5);
set(gcf,'name','Reconstructed components RCs')
clf;
for m=1:4
  subplot(4,2,2*m-1);
  plot(RC1(:,m),'r-');
  ylabel(sprintf('RC %d',m));
  ylim([-1 1]);

  subplot(4,2,2*m);
  plot(RC2(:,m),'r-');
  ylabel(sprintf('RC %d',m));
  ylim([-1 1]);
end

figure(6);
set(gcf,'name','Original time series X and reconstruction RC')
clf;
subplot(2,2,1)
plot(t,dop,'b-',t,sum(RC1,2),'r-');
subplot(2,2,2)
plot(t,vi,'b-',t,sum(RC2,2),'r-');
legend('Original','full reconstruction');

subplot(2,2,3)
plot(t,dop,'b',t,sum(RC1(:,1:2),2),'r');
subplot(2,2,4)
plot(t,vi,'b',t,sum(RC2(:,1:2),2),'r');
legend('Original','RCs 1-2');

%forecasting of RC1
%T represents the length of the time series
T = length(X) 
%A represents a zero matrix 
A = zeros(M-1,1);
for i = 1:rep %change
    A = A + RHO(M,i)*RHO(1:M-1,i);
end
v = norm(RHO(M,1:rep));%change 
%v represents the normalized form of eigen vectors 
A = A/(1-v^2);
A = flipdim(A,1);
%G represents the sum of the RC1 components 
G = sum(RC1(:,1:rep),2);%change
F = [G.', zeros(1,fop)].';
for i = T+1:T+fop
    for j = 1:M-1
        F(i) = F(i) + A(j)*F(i-j);
    end
end
F = F(T+1:T+fop);
%F represents the normalized form of the forecasted datas of dop
F = [G.', F.'].'; 
xlswrite('outputfile',F,'wbF');

%forecasting of RC2 components 
%B represents a zero matrix 
B = zeros(M-1,1);
for i = 1:rep %change
    B = B + RHO(M,i)*RHO(1:M-1,i);
end
v2 = norm(RHO(M,1:rep));%change
%v represents the normalized form of eigen vectors 
B = B/(1-v2^2);
B = flipdim(B,1);
I = sum(RC2(:,1:rep),2);%change
%I represents the sum of the RC2 components 
H = [I.', zeros(1,fop)].';
for i = T+1:T+fop
    for j = 1:M-1
        H(i) = H(i) + B(j)*H(i-j);
    end
end
H = H(T+1:T+fop);
H = [I.', H.'].' ;
%H represents the forecasted values of rh
xlswrite('outputfile',H,'VISF');

%denormalizing both the time series by multiplying F and H by standard
%deviation and adding mean of the respective time series 
dpd = F*stddop+meandop ;
rhd = H*stdrh+meanrh ;
xlswrite('outputfile',dpd,'wbFDN');
xlswrite('outputfile',rhd,'VISFDN');
