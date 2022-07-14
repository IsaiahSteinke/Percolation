C     ==================================================================
C     Program: bond_perc.f
C     Author: Isaiah Steinke
C
C     Runs multiple trials of the bond percolation code on a
C     triangular lattice in two dimensions.  Details for the operation
C     of the basic bond percolation code are given in the program
C     "bond.f".  In this code, the user specifies the lattice geometry
C     and the number of trials to run.  Each trial is run until a
C     lattice-spanning (or "infinite cluster") is present.
C
C     The code outputs a single text file of data containing the
C     following:
C        -the random number used to seed each trial (for debugging
C         purposes)
C        -the fraction of sites filled when each trial was terminated
C        -the number of bonds (i.e. the size) in the largest cluster
C        -the size of the infinite cluster.
C
C     *Critical percolation threshold for the bond problem on an
C     infinite triangular lattice is 0.3473.
C
C     Last Revision: July 8, 2010
C     ==================================================================

      program bond_perc

      implicit none
      integer m, n, t, pbc, rc, nb, bf, topfill, botfill
      integer scn, bcn, maxcn, maxcs, cln, lcn, lcs, rlc, clsum
      integer i, j, k, l, ii, seed, perccln, perccls, numtrials
      integer b(1000000,3), border(1000000,2), c(1000000)
      integer nn(10), btemp(2), nnb(25,4), tseed(50000)
      real rand
      double precision fb
      external nearestn
      common m, n, t, pbc, nn, scn

      open(unit=10, file='bond_perc.txt')

C     lattice dimensions, m x n

         m = 50                   !number of lattice sites wide
         n = 50                   !number of lattice sites high
         t = m*n                  !total no. of lattice sites

C     toggle for periodic boundary conditions on the left/right sides
C     note that the pbc flag only makes sense for m even with the
C     current lattice geometry (triangular)

         pbc = 0                  !set to 1 if using p.b.c.

C     lattice-specific geometry constants (need to be changed if the
C     lattice is changed)

         scn = 6                  !site coordination number
         bcn = 10                 !bond coordination number

C     calculate the overall total number of bonds possible (depends on
C     lattice geometry)

         if (pbc.eq.1) then
            nb = m*(3*n-2)
         else
            nb = (3*m*n)-(2*m)-(2*n)+1
         end if

C     number of percolation trials to execute and random number seed

         numtrials = 10           !upper limit: declared array size for tseed
         seed = 58302             !must be an integer
         call srand(seed)         !seed random number generator
C        generate the array of random number seeds for each trial
         do i=1,50000             !second number must be equal to declared array size
            tseed(i) = int(rand(0)*1000000)+1
         end do

C     completely zero out the nearest neighbor array

         do i=1,10
            nn(i) = 0             !nearest neighbor array
         end do

C     ==================================================================
C     main program loop
C     ==================================================================

         do ii=1,numtrials

C     initialize arrays to zero

         do j=1,3
            do i=1,1000000
               b(i,j) = 0         !bond array/matrix
            end do
         end do

         do j=1,2
            do i=1,1000000
               border(i,j) = 0    !bond order array
            end do
         end do

         do i=1,1000000
            c(i) = 0              !cluster array
         end do

C     initialize constants for checking for the infinite cluster

         maxcn = 0                !number of the largest overall cluster
         maxcs = 0                !size of the largest overall cluster

C     ------------------------------------------------------------------
C     generate the list of bonds in the bond array and set the
C     occupation equal to zero; initialize the bond order array
C     ------------------------------------------------------------------

         rc = 1                   !row number (counter)

         do i = 1, (t-1)
            call nearestn(i)
            do j = 1, scn
               if (nn(j).gt.i) then
                  b(rc,1) = i
                  b(rc,2) = nn(j)
                  b(rc,3) = 0
                  border(rc,1) = i
                  border(rc,2) = nn(j)
               else
                  goto 310
               end if
               rc = rc + 1
310            continue
            end do
         end do
C     ------------------------------------------------------------------

C     ------------------------------------------------------------------
C     randomize the order in which the bonds are filled in the bond
C     order array
C     ------------------------------------------------------------------

         call srand(tseed(ii))    !seeds the random number generator

C        permute bonds randomly in the bond order array

         do i = 1, nb
            j = i + (nb-i+1)*rand(0)
            btemp(1) = border(i,1)
            btemp(2) = border(i,2)
            border(i,1) = border(j,1)
            border(i,2) = border(j,2)
            border(j,1) = btemp(1)
            border(j,2) = btemp(2)
         end do
C     ------------------------------------------------------------------

C     ------------------------------------------------------------------
C     fill lattice with sites to percolation
C     ------------------------------------------------------------------

C     initialize constants/counters/calculated values

         bf = 0                   !number of bonds filled
         fb = real(bf)/real(nb)   !fraction of bonds filled
         cln = 1                  !lowest unused cluster number

         do i = 1,nb

C     initialize constants for tracking largest cluster characteristics

            lcn = 0               !cluster number of largest cluster
            lcs = 0               !size of largest cluster
            rlc = 0               !n.n. with largest cluster

C     initialize the nearest neighbor bond matrix to zero

            do j = 1,4
               do k = 1,25
                  nnb(k,j) = 0
               end do
            end do

C     ------------------------------------------------------------------
C     populate the nearest neighbor bond matrix with useful info
C     ------------------------------------------------------------------

            rc = 1                !reset row number/counter

            call nearestn(border(i,1))
            do j = 1,scn
               if (nn(j).ne.0) then
                  if (nn(j).ne.border(i,2)) then
                     if(nn(j).gt.border(i,1)) then
                        nnb(rc,1) = border(i,1)
                        nnb(rc,2) = nn(j)
                     else
                        nnb(rc,1) = nn(j)
                        nnb(rc,2) = border(i,1)
                     end if
                     rc = rc + 1
                  end if
               end if
            end do

            call nearestn(border(i,2))
            do j = 1,scn
               if (nn(j).ne.0) then
                  if (nn(j).ne.border(i,1)) then
                     if(nn(j).gt.border(i,2)) then
                        nnb(rc,1) = border(i,2)
                        nnb(rc,2) = nn(j)
                     else
                        nnb(rc,1) = nn(j)
                        nnb(rc,2) = border(i,2)
                     end if
                     rc = rc + 1
                  end if
               end if
            end do

            do j = 1,nb
               do k = 1,bcn
                  if(b(j,1).eq.nnb(k,1)) then
                     if(b(j,2).eq.nnb(k,2)) then
                        nnb(k,3) = b(j,3)
                        nnb(k,4) = c(b(j,3))
                     end if
                  end if
               end do
            end do
C     ------------------------------------------------------------------

C     find the nearest neighbor bond attached to the largest cluster
C     and get the largest cluster number and cluster size associated
C     with that nearest neighbor bond

            lcn = nnb(1,3)        !current largest cluster number
            lcs = nnb(1,4)        !current largest cluster size
            rlc = 1               !row of data in the n.n. bond matrix containing the largest cluster

            do k = 2,bcn
               if (nnb(k,1).ne.0) then           !no valid n.n. bond if =0
                  if (nnb(k,3).ne.0) then        !n.n. unoccupied if =0
                     if (nnb(k,4).gt.lcs) then
                        lcn = nnb(k,3)
                        lcs = nnb(k,4)           !or lcs = c(nnb(k,3))
                        rlc = k
                     end if
                  end if
               end if
            end do

C     Case 1: no nearest neighbors (i.e. lcs = 0)

            if (lcs.eq.0) then
               do j = 1,nb
                  if (b(j,1).eq.border(i,1)) then
                     if (b(j,2).eq.border(i,2)) then
                        b(j,3) = cln
                     end if
                  end if
               end do
               c(cln) = 1
               cln = cln+1
               goto 320
            end if

C     Case 2: one or more nearest neighbors occupied (i.e. lcs > 0)

            clsum = lcs                                    !temp variable for summing cluster size

            do k = 1,bcn
               if (nnb(k,1).ne.0) then                     !no valid n.n. if =0 (i.e. corner/side site)
                  if (nnb(k,3).ne.0) then                  !n.n. unoccupied if =0
                     if (nnb(k,3).ne.lcn) then             !filters out n.n. already belonging to the largest cluster
                        if (k.ne.1) then                   !this should check to make sure the same cluster number
                           do l = 1,(k-1)                  !is not "counted" more than once
                              if (nnb(l,3).eq.nnb(k,3)) then
                                 goto 330
                              end if
                           end do
                        end if
                        clsum = clsum + nnb(k,4)           !adds cluster size to largest cluster
                        do j = 1,nb
                           if (b(j,3).eq.nnb(k,3)) then    !this loop re-numbers all the bonds belonging to the smaller
                              b(j,3) = lcn                 !cluster with the largest cluster number
                           end if
                        end do
330                     c(nnb(k,3)) = 0                    !smaller cluster now contains zero sites
                     end if
                  end if
               end if
            end do

C     add randomly picked bond to the largest cluster

            do j = 1,nb
               if (b(j,1).eq.border(i,1)) then
                  if (b(j,2).eq.border(i,2)) then
                     b(j,3) = lcn
                  end if
               end if
            end do
            clsum = clsum+1
            c(lcn) = clsum

320         bf = bf + 1
            fb = real(bf)/real(nb)

C     this small loop should check for the largest overall cluster
C     and update the variables for the overall largest cluster size
C     and number

            if (c(lcn).gt.maxcs) then
               maxcs = c(lcn)
               maxcn = lcn
C           this sub-loop fixes the case where only one bond is occupied
            else if (lcs.eq.0) then
               if (maxcs.eq.0) then
                  maxcs = 1
                  maxcn = 1
               end if
            end if

C     check for percolation, i.e. presence of lattice spanning cluster

         if (i.ge.(n-1)) then     !need to occupy at least n-1 bonds before infinite cluster is possible
            do l = 1,(cln-1)
               botfill = 0        !used for checking if the bottom row is connected to the cluster
               topfill = 0        !used for checking if the top row is connected to the cluster
               if (c(l).ge.(n-1)) then   !checks for clusters with appropriate size that could span the lattice
                  do j = 1,m
                     do k = 1,nb
                        if (b(k,1).eq.j) then
                           if (b(k,3).eq.l) then
                              botfill = 1
                              goto 340
                           end if
                        end if
                     end do
                  end do
                  goto 360
340               do j = (t-m+1),t
                     do k = 1,nb
                        if (b(k,2).eq.j) then
                           if (b(k,3).eq.l) then
                              topfill = 1
                              goto 350
                           end if
                        end if
                     end do
                  end do
350               if ((topfill+botfill).eq.2) then
                     perccln = l     !infinite cluster number
                     perccls = c(l)  !infinite cluster size
                     goto 370
                  end if
360            end if
            end do
         end if

         end do

C     output information when infinite cluster is present
370      write(6,*) tseed(ii), fb, maxcs, perccls
         write(10,111) tseed(ii), fb, maxcs, perccls

111      format(i10,",",f12.9,",",i10,",",i10)

         end do

C     ==================================================================
C     end main program loop
C     ==================================================================

         stop
      end

C     ==================================================================
C     nearest neighbor subroutine that should populate the the nearest
C     neighbor array with the site numbers of the nearest neighbors for
C     a particular site - needs to be changed for different lattice
C     geometry
C     ==================================================================

      subroutine nearestn(rn)

      implicit none
      integer rn, m, n, t, pbc, scn
      integer z
      integer nn(10)
      common m, n, t, pbc, nn, scn

C     zero out the nearest neighbor array

         do z = 1,scn
            nn(z) = 0
         end do

C     check: lower-left corner site
            if (rn.eq.1) then
               nn(1) = rn+1
               nn(2) = rn+m
               nn(3) = rn+(m+1)
               if (pbc.eq.1) then
                  nn(4) = rn+(m-1)
                  nn(5) = rn+(2*m-1)
               end if
               goto 210
            end if

C     check: lower-right corner site
            if (rn.eq.m) then
               nn(1) = rn-1
               nn(2) = rn+m
               if (mod(m,2).eq.1) then !check for odd m
                  nn(3) = rn+(m-1)
                  goto 210
               end if
               if (pbc.eq.1) then
                  nn(3) = 1
               end if
               goto 210
            end if

C     check: upper-left corner site
            if (rn.eq.(t-(m-1))) then
               nn(1) = rn-m
               nn(2) = rn+1
               if (pbc.eq.1) then
                  nn(3) = t
               end if
               goto 210
            end if

C     check: upper-right corner site
            if (rn.eq.t) then
               if (mod(m,2).eq.1) then !check for odd m
                  nn(1) = rn-m
                  nn(2) = rn-1
                  goto 210
               end if
               nn(1) = rn-(m+1)
               nn(2) = rn-m
               nn(3) = rn-1
               if (pbc.eq.1) then
                  nn(4) = rn-(2*m-1)
                  nn(5) = rn-(m-1)
               end if
               goto 210
            end if

C     check: bottom row site (i.e. source edge)
            if (rn.lt.m) then
               if (mod(rn,2).eq.0) then
                  nn(1) = rn-1
                  nn(2) = rn+1
                  nn(3) = rn+m
               else
                  nn(1) = rn-1
                  nn(2) = rn+1
                  nn(3) = rn+(m-1)
                  nn(4) = rn+m
                  nn(5) = rn+(m+1)
               end if
               goto 210
            end if

C     check: top row site (i.e. drain edge)
            if (rn.gt.(t-m)) then
               if (mod(rn,2).eq.0) then
                  nn(1) = rn-(m+1)
                  nn(2) = rn-m
                  nn(3) = rn-(m-1)
                  nn(4) = rn-1
                  nn(5) = rn+1
               else
                  nn(1) = rn-m
                  nn(2) = rn-1
                  nn(3) = rn+1
               end if
               goto 210
            end if

C     check: left edge site
            if (mod((rn-1),m).eq.0) then
               nn(1) = rn-m
               nn(2) = rn+1
               nn(3) = rn+m
               nn(4) = rn+(m+1)
               if (pbc.eq.1) then
                  nn(5) = rn+(m-1)
                  nn(6) = rn+(2*m-1)
               end if
               goto 210
            end if

C     check: right edge site
            if (mod(rn,m).eq.0) then
               if (mod(m,2).eq.1) then !check for odd m
                  nn(1) = rn-m
                  nn(2) = rn-1
                  nn(3) = rn+(m-1)
                  nn(4) = rn+m
                  goto 210
               end if
               nn(1) = rn-(m+1)
               nn(2) = rn-m
               nn(3) = rn-1
               nn(4) = rn+m
               if (pbc.eq.1) then
                  nn(5) = rn-(2*m-1)
                  nn(6) = rn-(m-1)
               end if
               goto 210
            end if

C     general interior site
            if (mod(m,2).eq.1) then  !check for odd m
               if(mod((rn/m),2).eq.0) then  !even row
                  if (mod(rn,2).eq.0) then
                     nn(1) = rn-(m+1)
                     nn(2) = rn-m
                     nn(3) = rn-(m-1)
                     nn(4) = rn-1
                     nn(5) = rn+1
                     nn(6) = rn+m
                  else
                     nn(1) = rn-m
                     nn(2) = rn-1
                     nn(3) = rn+1
                     nn(4) = rn+(m-1)
                     nn(5) = rn+m
                     nn(6) = rn+(m+1)
                  end if
               else                         !odd row
                  if (mod(rn,2).eq.0) then
                     nn(1) = rn-m
                     nn(2) = rn-1
                     nn(3) = rn+1
                     nn(4) = rn+(m-1)
                     nn(5) = rn+m
                     nn(6) = rn+(m+1)
                  else
                     nn(1) = rn-(m+1)
                     nn(2) = rn-m
                     nn(3) = rn-(m-1)
                     nn(4) = rn-1
                     nn(5) = rn+1
                     nn(6) = rn+m
                  end if
               end if
            else
               if (mod(rn,2).eq.0) then
                  nn(1) = rn-(m+1)
                  nn(2) = rn-m
                  nn(3) = rn-(m-1)
                  nn(4) = rn-1
                  nn(5) = rn+1
                  nn(6) = rn+m
               else
                  nn(1) = rn-m
                  nn(2) = rn-1
                  nn(3) = rn+1
                  nn(4) = rn+(m-1)
                  nn(5) = rn+m
                  nn(6) = rn+(m+1)
               end if
            end if

210      end subroutine
C     ==================================================================
