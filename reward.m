function r_t=reward(t)
sigma=1;
r_t=1/2*exp( -(t-20).^2/(2*sigma^2) );
end