% bates and watts PCB example linear regression

function linearregression_weldingexample_withQR

figure(1); clf; figure(2);clf

x=[2 2.5 2.5 2.75 3 3 3]; x=x'; 
y=[89 97 91 98 100 104 97]; y=y';

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

Fvalue=finv(0.95,P,N-P);
xplot=1.8:0.1:3.2; xplot=xplot';
XPLOT=[ones(size(xplot)) xplot];
model=XPLOT*b;


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
params=[b-seB b+seB]

xlowlim=40; xhighlim=100;
ylowlim=0; yhighlim=25;
hold on;
plot([params(1,1) params(1,1)],[ylowlim yhighlim],'k--')
plot([params(1,2) params(1,2)],[ylowlim yhighlim],'k--')
plot([xlowlim xhighlim],[params(2,1) params(2,1)],'k--')
plot([xlowlim xhighlim],[params(2,2) params(2,2)],'k--')
plot([b(1)],b(2),'k+')

end
