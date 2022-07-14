% Matlab program used to plot output from site-bond percolation Fortran
% code.
% This code is written for the triangular lattice.
% Author: Isaiah Steinke
% Last Revision: November 18, 2010

% The widths of the lines and size of the points might need to be adjusted
% based upon lattice size for overall look.

% IMPORTANT: Output from the percolation code must be in the same directory
% as this Matlab program!

clear all;

% Specify lattice and program parameters.
m = 10;    % number of sites in the x-direction
n = 15;    % number of sites in the y-direction
a = 10;    % length of a bond
pbc = 0;   % periodic boundary condition flag
bs = 0;    % bondsite code flag, set to 1 if using the bondsite code instead of sitebond

% Parameter for distinguishing the lattice-spanning cluster (if different
% from the biggest cluster).  Be sure to set this variable to zero if you
% are not going to use this!  You may alternatively use this to highlight
% another cluster.
perccln = 0;   % cluster number of the lattice-spanning cluster (needs to be provided from code output)

% Calculate total number of bonds (dependent on lattice geometry).
if pbc == 0;
    numbonds = (3*m*n)-(2*m)-(2*n)+1;
else
    numbonds = m*(3*n-2);
end

% Specify bond and site coloring.
UnoccBondColor = [0.6 0.6 0.6];     % Unoccupied bond color (grey)
OccBondColor = [0.47 0.52 1.0];     % Occupied bond color, not belonging to the largest cluster (light blue)
OccSiteColor = [0.08 0.17 0.55];    % Occupied site color, not belonging to the largest cluster (dark blue)
MaxClBondColor = [0.25 0.75 0.25];  % Occupied bond color belonging to the largest cluster (green)
MaxClSiteColor = [0.11 0.31 0.21];  % Occupied site color belonging to the largest cluster (dark green)
PercBondColor = [0.75 0.25 0.25];   % Lattice-spanning bond color, if used (light red)
PercSiteColor = [0.6 0.0 0.0];      % Lattice-spanning site color, if used (red)

% Create (x,y) coordinate matrix for all of the sites (dependent on lattice
% geometry).
for i = 1:m   % first row of sites
    % check for site that is a multiple of m
    if mod(i,m) == 0
        S(i,1) = 0.5*sqrt(3)*a*m;  % x-coordinate
        if mod(i,2) == 0  % y-coordinate, need to check for even numbered site
            S(i,2) = 0.5*a;
        else
            S(i,2) = a;
        end
    else  % any other site that is not a multiple of m
        S(i,1) = 0.5*sqrt(3)*a*mod(i,m);  % x-coordinate
        if mod(i,2) == 0  % y-coordinate, need to check for even numbered site
            S(i,2) = 0.5*a;
        else
            S(i,2) = a;
        end
    end
end
for i = 1:(n-1)  % rest of the sites are just displaced multiples of a upwards
    for j = 1:m
        S((i*m+j),1) = S(j,1);        % x-coordinate
        S((i*m+j),2) = S(j,2)+(i*a);  % y-coordinate
    end
end

% Read in the list of all the bonds on the lattice and their occupation.
if (bs == 1)
    B = dlmread('bsbond.txt',',');
else
    B = dlmread('sbbond.txt',',');
end    

% Read in the site occupation and cluster data.
if (bs == 1)
    C = dlmread('bssite.txt',',');
else
    C = dlmread('sbsite.txt',',');
end    

% Get the largest cluster number.
if (bs == 1)
    for i = 1:numbonds
        if C(i,3) == max(C(:,3))
            lcn = i;  % largest cluster number
        end
    end
else
    for i = 1:(m*n)
        if C(i,3) == max(C(:,3))
            lcn = i;   % largest cluster number
        end
    end
end    

% Draw bonds connecting all the sites appropriately (does not draw periodic
% boundary bonds).
for i = 1:numbonds
    if ~((mod(B(i,1),m) == 1) & (mod(B(i,2),m) == 0))  % checks for non-pbc bond
        if (B(i,3) > 0)  % bond is occupied
            if (B(i,3) == perccln)  % bond is part of infinte cluster
                line([S(B(i,1),1) S(B(i,2),1)],[S(B(i,1),2) S(B(i,2),2)],'LineWidth',4,'Color',PercBondColor);
            elseif (B(i,3) == lcn)  % bond is part of largest cluster
                line([S(B(i,1),1) S(B(i,2),1)],[S(B(i,1),2) S(B(i,2),2)],'LineWidth',4,'Color',MaxClBondColor);
            else  % bond is not part of largest cluster
                line([S(B(i,1),1) S(B(i,2),1)],[S(B(i,1),2) S(B(i,2),2)],'LineWidth',4,'Color',OccBondColor);
            end
        else  % bond is unoccupied
            line([S(B(i,1),1) S(B(i,2),1)],[S(B(i,1),2) S(B(i,2),2)],'LineWidth',1,'Color',UnoccBondColor);
        end
    end
end
    
% Plot out the occupied sites.
hold on;
for i = 1:(m*n)
    if C(i,2) > 0  % checks for site occupation
        if (C(i,2) == perccln)  % site is part of infinte cluster
            plot(S(i,1),S(i,2),'ko','MarkerFaceColor',PercSiteColor,'MarkerEdgeColor','k','MarkerSize',6);
        elseif (C(i,2) == lcn)  % site is part of largest cluster
            plot(S(i,1),S(i,2),'ko','MarkerFaceColor',MaxClSiteColor,'MarkerEdgeColor','k','MarkerSize',6);
        else  % site is not part of largest cluster
            plot(S(i,1),S(i,2),'ko','MarkerFaceColor',OccSiteColor,'MarkerEdgeColor','k','MarkerSize',6);
        end
    end
end
axis equal;
axis([0 (0.5*sqrt(3)*a*(m+1)) 0 (a*(n+1))]);  % depends on lattice geometry
axis off;