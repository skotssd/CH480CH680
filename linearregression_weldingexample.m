% box hunter and hunter welding example linear regression

figure(1); clf; figure(2);clf

x=[2 2.5 2.5 2.75 3 3 3];
y=[89 97 91 98 100 104 97];
figure(1)
plot(x,y,'ko')

% linear regression
X=[x' ones(size(x'))]; y=y';
b=inv(X'*X)*X'*y

model=x*b(1)+b(2);

hold on

plot(x,model,'k')

% prediction interval

ssquared=sum((model-y').^2)./(size(x,2)-2); s=sqrt(ssquared);

for i=1:size(x,2)
	xh=[x(i); 1];
	upper(i)=xh'*b+s*sqrt(xh'*inv(X'*X)*xh)*sqrt(2*4.74);
	lower(i)=xh'*b-s*sqrt(xh'*inv(X'*X)*xh)*sqrt(2*4.74);
end
	
plot(x,upper,'k--');
plot(x,lower,'k--');


% parameter interval

Fvalue=finv(0.95,2,length(x)-2)
%Fvalue=4.74; 
contour95=2*ssquared*Fvalue;

intercept=40:0.4:90; slope=4:0.5:20;

for i=1:size(slope,2)
	for j=1:size(intercept,2)
		S(i,j)=sum((y'-(x*slope(i)+intercept(j))).^2);
	end
end

figure(2)
contour(intercept,slope,S,10); hold on;
contour(intercept,slope,S,[contour95 contour95])

% marginal estimates

invXtX=inv(X'*X);
se_beta1=s*sqrt(invXtX(1,1))
se_beta2=s*sqrt(invXtX(2,2))
tval=tinv(0.975,length(x)-2)
%tval=2.165;

figure(2); hold on

plot([40 90],[b(1) b(1)],'k--')
plot([b(2) b(2)],[4 20],'k--')
plot([b(2)+tval*se_beta2 b(2)+tval*se_beta2],[4 20],'k:')
plot([b(2)-tval*se_beta2 b(2)-tval*se_beta2],[4 20],'k:')
plot([40 90],[b(1)+tval*se_beta1 b(1)+tval*se_beta1],'k:')
plot([40 90],[b(1)-tval*se_beta1 b(1)-tval*se_beta1],'k:')


