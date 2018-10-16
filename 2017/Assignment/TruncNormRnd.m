function r = TruncNormRnd(mu,sigma,l,u)

% draw from a truncated normal distribution

% INPUT:
% mu: mean of the distribution
% sigma: standard deviation
% l: lower limit
% u: upper limit

% OUTPUT:
% r: value drawn from the truncated normal distribution 

t=truncate(makedist('normal',mu,sigma),l,u);

r=random(t,1);

end

