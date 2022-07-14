C     ==================================================================
C     Program: bondc.f
C     Author: Isaiah Steinke
C
C     Runs a bond percolation problem on a two dimensional triangular
C     lattice.  Bonds are randomly occupied from the empty lattice
C     state up to a specified fraction of the lattice, pb.  All
C     sites are assumed to be filled/occupied.  After filling the
C     lattice to the specified fraction of bonds, the presence of a
C     lattice-spanning (percolating) cluster is checked.  If a lattice-
C     spanning cluster is present, the conductance of the resulting
C     cluster is calculated.
C
C     This code primarily utilizes two arrays:
C        -the bond array (matrix): the two sites a bond connects is
C         listed in the first two columns and the cluster number the
C         bond belongs to is listed in the third column
C        -the cluster array: lists the size of the clusters in numerical
C         order
C
C     Another array is used to determine the order the bonds are
C     connected by randomly swapping bonds from an initial ordered list.
C
C     Code outputs three text files:
C        -bondorder.txt lists the order in which the bonds are
C         occupied
C        -bond.txt lists the contents of the bond array/matrix and the
C         cluster array
C        -bondocc.txt lists info for what happened when each bond was
C         occupied (useful for debugging).
C
C     *Critical percolation threshold for the bond problem on an
C     infinite triangular lattice is 0.3473.
C
C     This code uses subroutines to calculate the conductance from:
C     "Numerical Recipes in Fortran 77: The Art of Scientific
C     Computing", by William Press, et. al., Cambridge Press, 2001.
C     http://www.nr.com
C
C     Last Revision: October 18, 2010
C     ==================================================================

      program bondc

      implicit none
      integer m, n, t, pbc, rc, nb, bf, tbonds, topfill, botfill
      integer scn, bcn, maxcn, maxcs, cln, lcn, lcs, rlc, clsum
      integer i, j, k, l, seed, perccln, perccls
      integer b(1000000,3), border(1000000,2), c(1000000)
      integer nn(10), btemp(2), nnb(25,4), ija(20000), iter
      real rand
      double precision pb, fb
      double precision Va, g0, G(2500,2500), Itemp(2500), Iout(2500)
      double precision rowsum, Gtemp(2500,2500), Vint(2500)
      double precision V(2500), Itop, Ibot, Gtop, Gbot
      double precision sa(20000), err
      external nearestn, snrm
      common m, n, t, pbc, nn, scn
      COMMON /mat/ sa, ija

      open(unit=10, file='bond.txt')
      open(unit=11, file='bondocc.txt')
      open(unit=12, file='bondorder.txt')

C     lattice dimensions, m x n

         m = 50                   !number of lattice sites wide
         n = 50                   !number of lattice sites high
         t = m*n                  !total no. of lattice sites

C     toggle for periodic boundary conditions on the left/right sides
C     note that the pbc flag only makes sense for m even with the
C     current lattice geometry (triangular)

         pbc = 0                  !set to 1 if using p.b.c.

C     user-specified fraction of bonds to fill

         pb = 0.35d+00

C     Seed for the random number generator

         seed = 62703             !user defined seed, must be an integer

C     Conductance calculation parameters

         Va = 1.00d+00            !applied voltage [V]
         g0 = 1.00d+00            !base bond conductance [1/ohm]

C     lattice-specific geometry constants (need to be changed if the
C     lattice is changed)

         scn = 6                  !site coordination number
         bcn = 10                 !bond coordination number

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

         do i=1,10
            nn(i) = 0             !nearest neighbor array
         end do

         do i=1,1000000
            c(i) = 0              !cluster array
         end do

C     calculate the overall total number of bonds possible (depends on
C     lattice geometry)

         if (pbc.eq.1) then
            nb = m*(3*n-2)
         else
            nb = (3*m*n)-(2*m)-(2*n)+1
         end if

C     initialize constants for checking for the infinite cluster

         maxcn = 0                !number of the largest overall cluster
         maxcs = 0                !size of the largest overall cluster
         perccln = 0              !lattice-spanning cluster number
         perccls = 0              !lattice-spanning cluster size

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

         call srand(seed)         !seeds the random number generator

C     permute bonds randomly in the bond order array

         do i = 1, nb
            j = i + (nb-i+1)*rand(0)
            btemp(1) = border(i,1)
            btemp(2) = border(i,2)
            border(i,1) = border(j,1)
            border(i,2) = border(j,2)
            border(j,1) = btemp(1)
            border(j,2) = btemp(2)
         end do

C     output bond order array to file

         do i = 1, nb
            write(12,121) border(i,1), border(i,2)
         end do
C     ------------------------------------------------------------------

C     ==================================================================
C     main program loop
C     ==================================================================

C     initialize constants/counters/calculated values

         bf = 0                   !number of bonds filled
         fb = real(bf)/real(nb)   !fraction of bonds filled
         tbonds = pb*nb           !total number of bonds to fill
         cln = 1                  !lowest unused cluster number

         do i = 1,tbonds
            write(6,*) "bond chosen:", border(i,1), border(i,2)
            write(11,*) "bond chosen:", border(i,1), border(i,2)

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

C     output the contents of the n.n. bond matrix (for debugging)

C            write(6,*) "n.n. bond matrix"
C            write(11,*) "n.n. bond matrix"
C            do j = 1,bcn
C               write(6,*) (nnb(j,k), k = 1,4)
C               write(11,*) (nnb(j,k), k = 1,4)
C            end do
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

C            if (lcs.ne.0) then
C               write(6,*) "n.n. bond w/largest cluster:", nnb(rlc,1),
C     &         nnb(rlc,2)
C               write(6,*) "largest cluster number:", lcn
C               write(6,*) "largest neighbor cluster size:", lcs
C               write(11,*) "n.n. bond w/largest cluster:", nnb(rlc,1),
C     &         nnb(rlc,2)
C               write(11,*) "largest cluster number:", lcn
C               write(11,*) "largest neighbor cluster size:", lcs
C            end if

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
               write(6,*) "*no n.n. occupied*"
               write(6,*) "bond assigned to cluster number", cln
               write(11,*) "*no n.n. occupied*"
               write(11,*) "bond assigned to cluster number", cln
               cln = cln+1
               goto 320
            end if

C     Case 2: one or more nearest neighbors occupied (i.e. lcs > 0)

            clsum = lcs                                    !temp variable for summing cluster size
            write(6,*) "*one or more n.n. occupied*"
            write(11,*) "*one or more n.n. occupied*"

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
C                write(6,*) "adding", nnb(k,4), " to largest cluster"
C                write(6,*) "largest cluster is now", clsum
C                write(11,*) "adding", nnb(k,4), " to largest cluster"
C                write(11,*) "largest cluster is now", clsum
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
            write(6,*) "bond assigned to cluster number", lcn
            write(6,*) "size of cluster number", lcn, " is now", c(lcn)
            write(11,*) "bond assigned to cluster number", lcn
            write(11,*) "size of cluster number", lcn, " is now", c(lcn)

320         bf = bf + 1
            fb = real(bf)/real(nb)

            write(6,*) "fraction of lattice filled:", fb
            write(6,*) "--------------------"
            write(11,*) "fraction of lattice filled:", fb
            write(11,*) "--------------------"

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

         end do

C     ==================================================================
C     end main program loop
C     ==================================================================

         write(6,*)
         write(6,*) "******************************"
         write(11,*)
         write(11,*) "******************************"

C     output the largest cluster number and size

         write(6,*) "largest overall cluster number:", maxcn
         write(6,*) "largest overall cluster size:", maxcs
         write(11,*) "largest overall cluster number:", maxcn
         write(11,*) "largest overall cluster size:", maxcs

C     check for percolation, i.e. presence of lattice spanning cluster

         do i = 1,(cln-1)
            botfill = 0           !used for checking if the bottom row is connected to the cluster
            topfill = 0           !used for checking if the top row is connected to the cluster
            if (c(i).ge.(n-1)) then   !checks for clusters with appropriate size that could span the lattice
               write(6,*) "testing cluster:", i   !useful debugging message
               write(11,*) "testing cluster:", i
               do j = 1,m
                  do k = 1,nb
                     if (b(k,1).eq.j) then
                        if (b(k,3).eq.i) then
                           botfill = 1
                           goto 340
                        end if
                     end if
                  end do
               end do
               write(6,*) "source end not connected"
               write(11,*) "source end not connected"
               goto 360
340            do j = (t-m+1),t
                  do k = 1,nb
                     if (b(k,2).eq.j) then
                        if (b(k,3).eq.i) then
                           topfill = 1
                           goto 350
                        end if
                     end if
                  end do
               end do
               write(6,*) "drain end not connected"
               write(11,*) "drain end not connected"
350            if ((topfill+botfill).eq.2) then
                  write(6,*) "infinite cluster present"
                  write(11,*) "infinite cluster present"
                  perccln = i     !infinite cluster number
                  perccls = c(i)  !infinite cluster size
                  write(6,*) "infinite cluster number:", perccln
                  write(11,*) "infinite cluster number:", perccln
                  write(6,*) "infinite cluster size:", perccls
                  write(11,*) "infinite cluster size:", perccls
                  goto 370
               end if
360         end if
         end do

         write(6,*) "no infinite cluster present"
         write(11,*) "no infinite cluster present"

370      write(6,*) "******************************"
         write(11,*) "******************************"

C     calculate the conductance of the infinite cluster (if present)
         if (perccln.gt.0) then
            do i = 1,2500       !zero out conductance matrices
               do j = 1,2500
                  G(i,j) = 0.00d+00
                  Gtemp(i,j) = 0.00d+00
               end do
            end do
            do i = 1,2500       !zero out vectors used during calcuation
               Itemp(i) = 0.00d+00
               Vint(i) = 0.00d+00
               Iout(i) = 0.00d+00
               V(i) = 0.00d+00
            end do
            Itop = 0.00d+00
            Ibot = 0.00d+00
            Gtop = 0.00d+00
            Gbot = 0.00d+00
            do i = 1,nb           !populate off-diagonal elements of cond. matrix
               if (b(i,3).eq.perccln) then
                  G(b(i,1),b(i,2)) = -g0
                  G(b(i,2),b(i,1)) = G(b(i,1),b(i,2))
               else
                  G(b(i,1),b(i,2)) = -1.00d-12
                  G(b(i,2),b(i,1)) = G(b(i,1),b(i,2))
               end if
               if (b(i,1).gt.(t-2*m)) then !populate "temporary current" vector
                  if (b(i,1).le.(t-m)) then
                     if (b(i,2).gt.(t-m)) then
                        Itemp((b(i,1)-m)) = Itemp((b(i,1)-m))
     &                  -(G(b(i,1),b(i,2))*Va)
                     end if
                  end if
               end if
            end do
            do i = 1,t            !populate diagonal elements of cond. matrix
               rowsum = 0
               do j = 1,t
                  rowsum = rowsum + G(i,j)
               end do
               G(i,i) = -rowsum
            end do
C            write(6,*)'G'
C            write(11,*)'G'
C            do i = 1,t
C               write(6,*) (G(i,j), j = 1,t)
C               write(11,*) (G(i,j), j = 1,t)
C            end do
C            write(6,*) "--------------------"
C            write(11,*) "--------------------"
C            write(6,*) 'Itemp'
C            write(11,*) 'Itemp'
C            write(6,*) (Itemp(i), i = 1,(t-2*m))
C            write(11,*) (Itemp(i), i = 1,(t-2*m))
C            write(6,*) "--------------------"
C            write(11,*) "--------------------"
            do i = 1,(t-2*m)      !populate Gtemp matrix properly
               do j = 1,(t-2*m)
                  Gtemp(i,j) = G((i+m),(j+m))
               end do
            end do
C            write(6,*) 'Gtemp'
C            write(11,*) 'Gtemp'
C            do i = 1,(t-2*m)
C               write(6,*) (Gtemp(i,j), j = 1,(t-2*m))
C               write(11,*) (Gtemp(i,j), j = 1,(t-2*m))
C            end do
C            write(6,*) "--------------------"
C            write(11,*) "--------------------"
C           convert Gtemp matrix to row-indexed sparse storage mode
            do i = 1,20000
               sa(i) = 0.00d+00
               ija(i) = 0
            end do
            call sprsin(Gtemp,(t-2*m),2500,1.00d-16,20000,sa,ija)
            write(6,*) "Calculating internal node voltages"
            write(11,*) "Calculating internal node voltages"
C           set constants for linbcg routine
            iter = 0
            err = 1.00d-02
C           calculate internal node voltages using linbcg
            call linbcg((t-2*m),Itemp,Vint,2,1.00d-08,2500,iter,err)
C            write(6,*) "--------------------"
C            write(11,*) "--------------------"
C            write(6,*) 'Vint'
C            write(11,*) 'Vint'
C            write(6,*) (Vint(i), i = 1,(t-2*m))
C            write(11,*) (Vint(i), i = 1,(t-2*m))
C            write(6,*) "--------------------"
C            write(11,*) "--------------------"
            do i = 1,t            !populate voltage vector with the correct voltages
               if (i.le.m) then
                  V(i) = 0.00d+00
               else
                  if (i.gt.(t-m)) then
                     V(i) = Va
                  else
                     V(i) = Vint((i-m))  !Vint contains the solutions from linbcg routine
                  end if
               end if
            end do
C            write(6,*) 'V'
C            write(11,*) 'V'
C            write(6,*) (V(i), i = 1,t)
C            write(11,*) (V(i), i = 1,t)
C            write(6,*) "--------------------"
C            write(11,*) "--------------------"
C           convert G matrix to row-indexed sparse storage mode
            do i = 1,20000
               sa(i) = 0.00d+00
               ija(i) = 0
            end do
            call sprsin(G,t,2500,1.00d-10,20000,sa,ija)
            write(6,*) "Calculating currents"
            write(11,*) "Calculating currents"
C           calculate currents by multiplying G*V
            call dsprsax(sa,ija,V,Iout,t)
C            write(6,*) 'Iout'
C            write(11,*) 'Iout'
C            write(6,*) (Iout(i), i = 1,t)
C            write(11,*) (Iout(i), i = 1,t)
            write(6,*) "--------------------"
            write(11,*) "--------------------"
            do i = 1,m
               Ibot = Ibot + Iout(i)  !Iout contains solution from dsprsax routine
               Itop = Itop + Iout((i+t-m))
            end do
            Gtop = Itop/Va
            Gbot = abs(Ibot)/Va
            write(6,*) "Conductance:", Gtop, Gbot
            write(11,*) "Conductance:", Gtop, Gbot
         end if


C     output bond and cluster arrays

         do j = 1,nb
            write(10,111) b(j,1), b(j,2), b(j,3), j, c(j)
         end do

111      format(i10,",",i10,",",i10,",",i10,",",i10)
121      format(i10,",",i10)

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

C     ==================================================================
C     Subroutines taken from Numerical Recipes...
C     ==================================================================

C     sprsin: converts a matrix to sparse storage
      SUBROUTINE sprsin(a,n,np,thresh,nmax,sa,ija)
      INTEGER n,nmax,np,ija(nmax)
      DOUBLE PRECISION thresh,a(np,np),sa(nmax)
      INTEGER i,j,k
      do 11 j=1,n
        sa(j)=a(j,j)
11    continue
      ija(1)=n+2
      k=n+1
      do 13 i=1,n
        do 12 j=1,n
          if(abs(a(i,j)).ge.thresh)then
            if(i.ne.j)then
              k=k+1
              if(k.gt.nmax)pause 'nmax too small in sprsin'
              sa(k)=a(i,j)
              ija(k)=j
            endif
          endif
12      continue
        ija(i+1)=k+1
13    continue
      return
      END
      
C     linbcg: solves a sparse matrix system of linear equations Ax=b
C     using the preconditioned biconjugate gradient method
      SUBROUTINE linbcg(n,b,x,itol,tol,itmax,iter,err)
      INTEGER iter,itmax,itol,n,NMAX
      DOUBLE PRECISION err,tol,b(*),x(*),EPS
      PARAMETER (NMAX=20000,EPS=1.00d-14)
CU    USES atimes,asolve,snrm
      INTEGER j
      DOUBLE PRECISION ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,
     *znrm,p(NMAX),pp(NMAX),r(NMAX),rr(NMAX),z(NMAX),zz(NMAX),snrm
      iter=0
      call atimes(n,x,r,0)
      do 11 j=1,n
        r(j)=b(j)-r(j)
        rr(j)=r(j)
11    continue
C     call atimes(n,r,rr,0)
      znrm=1.00d+00
      if(itol.eq.1) then
        bnrm=snrm(n,b,itol)
      else if (itol.eq.2) then
        call asolve(n,b,z,0)
        bnrm=snrm(n,z,itol)
      else if (itol.eq.3.or.itol.eq.4) then
        call asolve(n,b,z,0)
        bnrm=snrm(n,z,itol)
        call asolve(n,r,z,0)
        znrm=snrm(n,z,itol)
      else
        pause 'illegal itol in linbcg'
      endif
      call asolve(n,r,z,0)
100   if (iter.le.itmax) then
        iter=iter+1
        zm1nrm=znrm
        call asolve(n,rr,zz,1)
        bknum=0.00d+00
        do 12 j=1,n
          bknum=bknum+z(j)*rr(j)
12      continue
        if(iter.eq.1) then
          do 13 j=1,n
            p(j)=z(j)
            pp(j)=zz(j)
13        continue
        else
          bk=bknum/bkden
          do 14 j=1,n
            p(j)=bk*p(j)+z(j)
            pp(j)=bk*pp(j)+zz(j)
14        continue
        endif
        bkden=bknum
        call atimes(n,p,z,0)
        akden=0.00d+00
        do 15 j=1,n
          akden=akden+z(j)*pp(j)
15      continue
        ak=bknum/akden
        call atimes(n,pp,zz,1)
        do 16 j=1,n
          x(j)=x(j)+ak*p(j)
          r(j)=r(j)-ak*z(j)
          rr(j)=rr(j)-ak*zz(j)
16      continue
        call asolve(n,r,z,0)
        if(itol.eq.1.or.itol.eq.2)then
          znrm=1.00d+00
          err=snrm(n,r,itol)/bnrm
        else if(itol.eq.3.or.itol.eq.4)then
          znrm=snrm(n,z,itol)
          if(abs(zm1nrm-znrm).gt.EPS*znrm) then
            dxnrm=abs(ak)*snrm(n,p,itol)
            err=znrm/abs(zm1nrm-znrm)*dxnrm
          else
            err=znrm/bnrm
            goto 100
          endif
          xnrm=snrm(n,x,itol)
          if(err.le.0.50d+00*xnrm) then
            err=err/xnrm
          else
            err=znrm/bnrm
            goto 100
          endif
        endif
        write (*,*) ' iter=',iter,' err=',err
      if(err.gt.tol) goto 100
      endif
      return
      END
      
C     atimes
      SUBROUTINE atimes(n,x,r,itrnsp)
      INTEGER n,itrnsp,ija,NMAX
      DOUBLE PRECISION x(n),r(n),sa
      PARAMETER (NMAX=20000)
      COMMON /mat/ sa(NMAX),ija(NMAX)
      if (itrnsp.eq.0) then
        call dsprsax(sa,ija,x,r,n)
      else
        call dsprstx(sa,ija,x,r,n)
      endif
      return
      END
      
C     asolve
      SUBROUTINE asolve(n,b,x,itrnsp)
      INTEGER n,itrnsp,ija,NMAX,i
      DOUBLE PRECISION x(n),b(n),sa
      PARAMETER (NMAX=20000)
      COMMON /mat/ sa(NMAX),ija(NMAX)
      do 11 i=1,n
        x(i)=b(i)/sa(i)
11    continue
      return
      END
      
C     snrm
      FUNCTION snrm(n,sx,itol)
      INTEGER n,itol,i,isamax
      DOUBLE PRECISION sx(n),snrm
      if (itol.le.3)then
        snrm=0.00d+00
        do 11 i=1,n
          snrm=snrm+sx(i)**2
11      continue
        snrm=sqrt(snrm)
      else
        isamax=1
        do 12 i=1,n
          if(abs(sx(i)).gt.abs(sx(isamax))) isamax=i
12      continue
        snrm=abs(sx(isamax))
      endif
      return
      END
      
C     dsprsax: multiplies a sparse matrix times a vector, A*x
      SUBROUTINE dsprsax(sa,ija,x,b,n)
      INTEGER n,ija(*)
      DOUBLE PRECISION b(n),sa(*),x(n)
      INTEGER i,k
      if (ija(1).ne.n+2) pause 'mismatched vector and matrix in dsprsax'
      do 12 i=1,n
        b(i)=sa(i)*x(i)
        do 11 k=ija(i),ija(i+1)-1
          b(i)=b(i)+sa(k)*x(ija(k))
11      continue
12    continue
      return
      END
      
C     dsprstx: multiplies the transpose of a matrix times a vector, A'*x
      SUBROUTINE dsprstx(sa,ija,x,b,n)
      INTEGER n,ija(*)
      DOUBLE PRECISION b(n),sa(*),x(n)
      INTEGER i,j,k
      if (ija(1).ne.n+2) pause 'mismatched vector and matrix in dsprstx'
      do 11 i=1,n
        b(i)=sa(i)*x(i)
11    continue
      do 13 i=1,n
        do 12 k=ija(i),ija(i+1)-1
          j=ija(k)
          b(j)=b(j)+sa(k)*x(i)
12      continue
13    continue
      return
      END
      
