function K = K(a, th)
%K(a,th) calculates the K point of a graphene lattice
%   Detailed explanation goes here
sinth = sin(th); 
costh = cos(th);
K = [costh, -sinth; sinth, costh] * [2*pi/3/a; 2*pi/3/sqrt(3)/a];

end

