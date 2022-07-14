C     ==================================================================
C     Program: bond.f
C     Author: Isaiah Steinke
C
C     Runs a bond percolation problem on a two dimensional square
C     lattice.  Bonds are randomly occupied from the empty lattice
C     state up to a specified fraction of the lattice, pb.  All
C     sites are assumed to be filled/occupied.
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
C     infinite square lattice is 1/2 (exact).
C
C     Last Revision: July 8, 2010
C     ==================================================================

      program bond
      
      implicit none
      integer m, n, t, pbc, rc, nb, bf, tbonds, topfill, botfill
      integer scn, bcn, maxcn, maxcs, cln, lcn, lcs, rlc, clsum
      integer i, j, k, l, seed, perccln, perccls
      integer b(1000000,3), border(1000000,2), c(1000000)
      integer nn(10), btemp(2), nnb(25,4)
      real rand
      double precision pb, fb
      external nearestn
      common m, n, t, pbc, nn, scn
      
      open(unit=10, file='bond.txt')
      open(unit=11, file='bondocc.txt')
      open(unit=12, file='bondorder.txt')

C     lattice dimensions, m x n

         m = 50                   !number of lattice sites wide
         n = 50                   !number of lattice sites high
         t = m*n                  !total no. of lattice sites
         
C     toggle for periodic boundary conditions on the left/right sides

         pbc = 0                  !set to 1 if using p.b.c.

C     user-specified fraction of bonds to fill

         pb = 0.50d+00

C     lattice-specific geometry constants (need to be changed if the
C     lattice is changed)

         scn = 4                  !site coordination number
         bcn = 6                  !bond coordination number

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
            nb = m*(2*n-1)
         else
            nb = (2*m*n)-m-n
         end if

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

         seed = 184489            !user defined seed, must be an integer
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

            write(6,*) "n.n. bond matrix"
            write(11,*) "n.n. bond matrix"
            do j = 1,bcn
               write(6,*) (nnb(j,k), k = 1,4)
               write(11,*) (nnb(j,k), k = 1,4)
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

            if (lcs.ne.0) then
               write(6,*) "n.n. bond w/largest cluster:", nnb(rlc,1),
     &         nnb(rlc,2)
               write(6,*) "largest cluster number:", lcn
               write(6,*) "largest neighbor cluster size:", lcs
               write(11,*) "n.n. bond w/largest cluster:", nnb(rlc,1),
     &         nnb(rlc,2)
               write(11,*) "largest cluster number:", lcn
               write(11,*) "largest neighbor cluster size:", lcs
            end if

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
               write(6,*) "site assigned to cluster number", cln
               write(11,*) "*no n.n. occupied*"
               write(11,*) "site assigned to cluster number", cln
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
                write(6,*) "adding", nnb(k,4), " to largest cluster"
                write(6,*) "largest cluster is now", clsum
                write(11,*) "adding", nnb(k,4), " to largest cluster"
                write(11,*) "largest cluster is now", clsum
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
            write(6,*) "site assigned to cluster number", lcn
            write(6,*) "size of cluster number", lcn, " is now", c(lcn)
            write(11,*) "site assigned to cluster number", lcn
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
               write(6,*) "testing cluster"   !useful debugging message
               write(11,*) "testing cluster"
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
               if (pbc.eq.1) then
                  nn(3) = m
               end if
               goto 210
            end if

C     check: lower-right corner site
            if (rn.eq.m) then
               nn(1) = rn-1
               nn(2) = rn+m
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
               nn(1) = rn-m
               nn(2) = rn-1
               if (pbc.eq.1) then
                  nn(3) = rn-(m-1)
               end if
               goto 210
            end if

C     check: bottom row site (i.e. source edge)
            if (rn.lt.m) then
               nn(1) = rn-1
               nn(2) = rn+1
               nn(3) = rn+m
               goto 210
            end if

C     check: top row site (i.e. drain edge)
            if (rn.gt.(t-m)) then
               nn(1) = rn-m
               nn(2) = rn-1
               nn(3) = rn+1
               goto 210
            end if

C     check: left edge site
            if (mod((rn-1),m).eq.0) then
               nn(1) = rn-m
               nn(2) = rn+1
               nn(3) = rn+m
               if (pbc.eq.1) then
                  nn(4) = rn+(m-1)
               end if
               goto 210
            end if

C     check: right edge site
            if (mod(rn,m).eq.0) then
               nn(1) = rn-m
               nn(2) = rn-1
               nn(3) = rn+m
               if (pbc.eq.1) then
                  nn(4) = rn-(m-1)
               end if
               goto 210
            end if

C     general interior site
            nn(1) = rn-m
            nn(2) = rn-1
            nn(3) = rn+1
            nn(4) = rn+m

210      end subroutine
C     ==================================================================
