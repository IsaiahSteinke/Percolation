% Matlab program used to plot output from bond percolation Fortran code.
% This code is written for the square lattice.
% Author: Isaiah Steinke
% Last Revision: August 11, 2010

% The widths of the lines and size of the points might need to be adjusted
% based upon lattice size for overall look.

% IMPORTANT: Output from the percolation code must be in the same directory
% as this Matlab program!

clear all;

% Specify lattice and program parameters.
m = 50;    % number of sites in the x-direction
n = 50;    % number of sites in the y-direction
a = 10;    % length of a bond
pbc = 0;   % periodic boundary condition flag
plotallsites = 0;  % set to 1 if you want all the sites plotted

% Parameter for distinguishing the lattice-spanning cluster (if different
% from the biggest cluster).  Be sure to set this variable to zero if you
% are not going to use this!  You may alternatively use this to highlight
% another cluster.
perccln = 0;  % cluster number of the infinite cluster (needs to be provided from code output)

% Calculate total number of bonds (dependent on lattice geometry).
if pbc == 0;
    numbonds = (2*m*n)-m-n;
else
    numbonds = m*(2*n-1);
end

% Specify bond coloring.
UnoccBondColor = [0.6 0.6 0.6];  % Unoccupied bond color (grey)
OccBondColor = [0.0 0.0 0.8];    % Occupied bond color, not belonging to largest cluster (blue)
MaxClBondColor = [0.0 0.8 0.0];  % Occupied bond color belonging to largest cluster (green)
PercBondColor = [0.8 0.0 0.0];   % Lattice-spanning bond color, if used (red)

% Create (x,y) coordinate matrix for all of the sites (dependent on lattice
% geometry).
for i = 1:m   % first row of sites
    % check for site that is a multiple of m
    if mod(i,m) == 0
        S(i,1) = a*m;         % x-coordinate
        S(i,2) = a;           % y-coordinate
    else  % any other site that is not a multiple of m
        S(i,1) = a*mod(i,m);  % x-coordinate
        S(i,2) = a;           % y-coordinate
    end
end
for i = 1:(n-1)  % rest of the sites are just displaced multiples of a upwards
    for j = 1:m
        S((i*m+j),1) = S(j,1);        % x-coordinate
        S((i*m+j),2) = S(j,2)+(i*a);  % y-coordinate
    end
end

% Read in bond occupation data from bond percolation code.
B = dlmread('bond.txt',',');

% Get the largest cluster number.
for i = 1:numbonds
    if B(i,5) == max(B(:,5))
        lcn = i;   % largest cluster number
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

% Plot out the sites.
hold on;
if plotallsites == 1
    plot(S(:,1),S(:,2),'ko','MarkerFaceColor','k','MarkerSize',6);
else
    for i = 1:numbonds
        if B(i,3) > 0
            plot(S(B(i,1),1),S(B(i,1),2),'ko','MarkerFaceColor','k','MarkerSize',6);
            plot(S(B(i,2),1),S(B(i,2),2),'ko','MarkerFaceColor','k','MarkerSize',6);
        end
    end
end
axis equal;
axis([0 (a*(m+1)) 0 (a*(n+1))]);  % depends on lattice geometry
axis off;