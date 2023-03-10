function PathLengths = pathlengths(lambda)
% in cm
Hb0 = 60*10^-6;
Hb = 40*10^-6;
g = 0.9;
c = 3*10^10;
e = GetExtinctions(lambda);
mua = e(:,1).*Hb0+e(:,2).*Hb;
mus = 150*(lambda/560).^(-2);
z0 = 1./((1-g)*mus');
gamma = sqrt(c./(3*(mua+(1-g)*mus')));
PathLengths = (c*z0./(2*gamma.*(mua*c).^0.5)).*(1+(3/c)*mua.*(gamma).^2);