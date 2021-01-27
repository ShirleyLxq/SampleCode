function H=bernoulliEntropy(p)
% a bernoulli distribution 
H=-p.*log2(p)-(1-p).*log2(1-p);