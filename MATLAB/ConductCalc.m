% Matlab program used to calculate the conductance of the lattice-spanning
% cluster.  The code will only consider conductances that belong to the
% infinte cluster, other bonds will have their conductances set to near 
% zero.  A potential is applied to the top row of sites and the bottom row 
% of sites are grounded.  The program will first calculate the potentials 
% at all the internal nodes and then calculate the total current at the top
% and bottom contacts.  Using this information, the conductance of the
% infinite cluster is calculated.
% The code has been tested and appears to function properly for all
% combinations of lattice and percolation types.
% Conductances may be universally fixed at some value or randomly assigned.
% Author: Isaiah Steinke
% Last Revision: August 24, 2010

% Note that the output from the percolation code should have an infinite
% cluster present!  It will either return bad values or (more likely) zero
% for the conductance.

% IMPORTANT: Output from the percolation code must be in the same directory
% as this Matlab program!

clear all;

% Specify lattice and parameters from Fortran code.
m = 10;       % lattice width (in # sites)
n = 10;       % lattice height (in # sites)
t = m*n;      % total number of lattice sites
pbc = 0;      % periodic boundary condition flag, set to 1 if used
perctype = 1; % percolation code type: 1:site, 2:bond, 3:site-bond
lattype = 1;  % lattice type: 1:square, 2:triangular
perccln = 1;  % cluster number of the lattice-spanning cluster (must be non-zero)

% Specify applied voltage conductance parameters.  The condtype flag will
% specify how the bond conductances are determined.  They can all be fixed
% to one value, g0, or be randomly distributed.
Va = 1;       % voltage applied to the top contact, units: [V]
g0 = 1;       % base conductance value for all the bonds, units: [1/ohm]
condtype = 1; % conductance type: 1:fixed, 2:variable
% The variable conductance is currently generated using a uniform
% distribution on [0,1) and multiplying this number by the base conductance
% value.

% Seed the random number generator if variable conductance is used.  The
% second argument (the number) seeds the random number generator.
if condtype == 2
    rand('twister',1838534);
end

% Calculate the total number of bonds.
if lattype == 1
    if pbc == 0
        numbonds = (2*m*n)-m-n;
    else
        numbonds = m*(2*n-1);
    end
end

if lattype == 2
    if pbc == 0
        numbonds = (3*m*n)-(2*m)-(2*n)+1;
    else
        numbonds = m*(3*n-2);
    end
end

% Read in the appropriate data from the percolation code output.
if perctype == 1
    B = dlmread('bondlist.txt',',');   % complete list of bonds
    C = dlmread('site.txt',',');       % site occupation and cluster data
end

if perctype == 2
    B = dlmread('bond.txt',',');       % bond list w/occupation and cluster data
end

if perctype == 3
    B = dlmread('sbbond.txt',',');     % list of bonds and their occupation
    C = dlmread('sbsite.txt',',');     % site occupation and cluster data
end

% Zero out conductance matrix and voltage and current vectors.
G = zeros(t,t);                        % conductance matrix
V = zeros(t,1);                        % vector of voltages
Itemp = zeros((t-2*m),1);              % vector used to calculate internal node voltages

% Use the percolation code data to populate the conductance matrix and the
% vector used to calculate the internal node voltages.
if perctype == 1
    for i = 1:numbonds
        if (C(B(i,1),2) == perccln) && (C(B(i,2),2) == perccln) % checks for bonds belonging to the infinite cluster
            if condtype == 1           % populate off-diagonal elements  
                G(B(i,1),B(i,2)) = -g0;
                G(B(i,2),B(i,1)) = G(B(i,1),B(i,2));
            else if condtype == 2
                    G(B(i,1),B(i,2)) = -g0*rand;
                    G(B(i,2),B(i,1)) = G(B(i,1),B(i,2));
                end
            end
        else
            G(B(i,1),B(i,2)) = -1E-12; % bonds not belonging to the infinite cluster have very low conductance
            G(B(i,2),B(i,1)) = -1E-12; % this helps to avoid singular matrices during inversion
        end
        if ((B(i,1) > (t-2*m)) && (B(i,1) <= (t-m))) % checks for nodes near the top row of sites
           if (B(i,2) > (t-m))
                Itemp(B(i,1)-m) = Itemp(B(i,1)-m)+(-G(B(i,1),B(i,2))*Va); % populate Itemp vector
            end
        end        
    end
end

if perctype == 2
    for i = 1:numbonds
        if (B(i,3) == perccln)         % checks for bonds belonging to the infinite cluster
            if condtype == 1           % populate off-diagonal elements  
                G(B(i,1),B(i,2)) = -g0;
                G(B(i,2),B(i,1)) = G(B(i,1),B(i,2));
            else if condtype == 2
                    G(B(i,1),B(i,2)) = -g0*rand;
                    G(B(i,2),B(i,1)) = G(B(i,1),B(i,2));
                end
            end
        else 
            G(B(i,1),B(i,2)) = -1E-12; % bonds not belonging to the infinite cluster have very low conductance
            G(B(i,2),B(i,1)) = -1E-12; % this helps to avoid singular matrices during inversion
        end
        if ((B(i,1) > (t-2*m)) && (B(i,1) <= (t-m))) % checks for nodes near the top row of sites
            if (B(i,2) > (t-m))
                Itemp(B(i,1)-m) = Itemp(B(i,1)-m)+(-G(B(i,1),B(i,2))*Va); % populate Itemp vector
            end
        end        
    end
end

if perctype == 3
    for i = 1:numbonds
        if (B(i,3) == perccln)         % checks for bonds belonging to the infinite cluster
            if (C(B(i,1),2) == perccln) && (C(B(i,2),2) == perccln) % checks if both end sites belong to inf. cluster
                if condtype == 1           % populate off-diagonal elements  
                    G(B(i,1),B(i,2)) = -g0;
                    G(B(i,2),B(i,1)) = G(B(i,1),B(i,2));
                else if condtype == 2
                        G(B(i,1),B(i,2)) = -g0*rand;
                        G(B(i,2),B(i,1)) = G(B(i,1),B(i,2));
                    end
                end
            else
                G(B(i,1),B(i,2)) = -1E-12; % dead end in the infinite cluster; make low conductance
                G(B(i,2),B(i,1)) = -1E-12; % this helps to avoid singular matrices during inversion
            end
        else
            G(B(i,1),B(i,2)) = -1E-12; % bonds not belonging to the infinite cluster have very low conductance
            G(B(i,2),B(i,1)) = -1E-12; % this helps to avoid singular matrices during inversion
        end
        if ((B(i,1) > (t-2*m)) && (B(i,1) <= (t-m))) % checks for nodes near the top row of sites
            if (B(i,2) > (t-m))
                Itemp(B(i,1)-m) = Itemp(B(i,1)-m)+(-G(B(i,1),B(i,2))*Va); % populate Itemp vector
            end
        end        
    end
end

rowsum = sum(G,2);      % returns a vector with the sum of all the elements in each row
for i = 1:t             % populate the diagonal elements
    G(i,i) = -(rowsum(i));
end

% Calculate the conductance.
Gtemp = G;              % copy conductance matrix into a new one
for i = 1:m             % delete the last m columns and rows
    Gtemp(:,(t-i+1)) = [];
    Gtemp((t-i+1),:) = [];
end
for i = 1:m             % delete the first m columns and rows
    Gtemp(:,1) = [];
    Gtemp(1,:) = [];
end
Vint = Gtemp\Itemp;     % solve for internal node voltages, supposedly better than inv(Gtemp)*Itemp
for i = 1:t             % rewrite voltage vector with the correct voltages
    if (i <= m)
        V(i) = 0;
    else if (i > (t-m))
            V(i) = Va;
        else
           V(i) = Vint(i-m);
        end
    end
end
I = G*V;                % calculate current vector
Ibot = 0;               % temp variable for summing currents (bottom)
Itop = 0;               % temp variable for summing currents (top)
for i = 1:m             % sum the currents going in/coming out of the bottom/top
    Ibot = I(i)+Ibot;
    Itop = I(t-i+1)+Itop;
end
Gbot = abs(Ibot)/Va;    % calculate the overall conductance using the total current into the bottom
Gtop = Itop/Va;         % calculate the overall conductance using the total current from the top

% ** Gbot should equal Gtop if everything is done correctly.**