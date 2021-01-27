function dopamineActivity=DA(x)
% x is the learning error delta(t), which can be a matrix or vector
xStar=0.27;
alpha=6; beta=6;

dopamineActivity=zeros(size(x));
isType1=x<0;
dopamineActivity(isType1)=x(isType1)/alpha;
isType2=and(x>=0, x<xStar);
dopamineActivity(isType2)=x(isType2);
isType3=x>=xStar;
dopamineActivity(isType3)=xStar+(x(isType3)-xStar)/beta;