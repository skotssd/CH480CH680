% bates and watts PCB example linear regression

function exercise4_L2_linearmodel_PCBdata
  
% for octave

pkg load statistics

figure(1); clf; figure(2);clf

[x,y]=getdata; x=x.^(1/3); y=log(y);

figure(1)
plot(x,y,'ko')

% linear regression
X=[ones(size(x)) x]; %y=y';
b=inv(X'*X)*X'*y;
P=2; % two parameters
N=length(x); %number of observatrions
nu=N-P; % degrees of freedom.  no. obs-no.parameters
bestmodel=X*b;
s2=sum((bestmodel-y).^2)./nu; s=sqrt(s2);

hold on

% 1-alpha condidence band for the response function for any x ------

Fvalue=finv(0.95,2,length(x)-2); xplot=0.9:0.1:2.4; 
XPLOT=[ones(size(xplot')) xplot']; model=XPLOT*b;

plot(xplot,model,'k')

for i=1:length(xplot)
	xh=[1; xplot(i)];
	upper(i)=xh'*b+s*sqrt(xh'*inv(X'*X)*xh)*sqrt(P*Fvalue);
	lower(i)=xh'*b-s*sqrt(xh'*inv(X'*X)*xh)*sqrt(P*Fvalue);
end

plot(xplot,upper,'k--')
plot(xplot,lower,'k--')

% contour of joint confidence interval. the ellipsoid

w=0:0.1:2.5*pi; % for plotting the elippse

[Q,R]=qr(X); R1 = R(1:2,1:2); invR1=inv(R1);

for i=1:length(w)
	    scalar=sqrt(P*s2*Fvalue);
        beta(:,i)=b+scalar*invR1*([cos(w(i)); sin(w(i))]);
    end

figure(2); plot(beta(1,:),beta(2,:))

% confidence intervals on parameter estimates

tvalue=tinv(0.975,length(x)-P);
seB=s*sqrt(diag(inv(X'*X)))*tvalue;
params=[b+seB b-seB];

xlowlim=-3.4; xhighlim=-1.4;
ylowlim=1.6; yhighlim=3;
hold on;
plot([params(1,1) params(1,1)],[ylowlim yhighlim],'k--')
plot([params(1,2) params(1,2)],[ylowlim yhighlim],'k--')
plot([xlowlim xhighlim],[params(2,1) params(2,1)],'k--')
plot([xlowlim xhighlim],[params(2,2) params(2,2)],'k--')
plot([b(1)],b(2),'k+')

end

function [x,y]=getdata

data=[...
    1   0.6
    1   1.6
    1   0.5
    1   1.2
    2   2.0
    2   1.3
    2   2.5
    3   2.2
    3   2.4
    3   1.2
    4   3.5
    4   4.1
    4   5.1
    5   5.7
    6   3.4
    6   9.7
    6   8.6
    7   4.0
    7   5.5
    7   10.5
    8   17.5
    8   13.4
    8   4.5
    9   30.4
    11  12.4
    12  13.4
    12  26.2
    12  7.4
];

x=data(:,1); y=data(:,2);

end
