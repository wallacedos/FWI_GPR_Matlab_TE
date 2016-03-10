function p = ricker( fr, t, center )
%RICKER Summary of this function goes here
%   Detailed explanation goes here

t = t - 1/fr * center;
p = (1 - 2*pi^2 * fr^2 .* t.^2) .* exp(- pi^2 * fr^2 .* t.^2);

% figure()
% plot(t,p)

end

