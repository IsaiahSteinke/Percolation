C     ==================================================================
C     Program: bs_perc.f
C     Author: Isaiah Steinke
C
C     Runs multiple site-bond problems on a square lattice.  The
C     user must specify an array of values for pb at which to run the
C     site-bond problem.  The code will occupy the specified fraction of
C     bonds and then occupy sites until a lattice-spanning cluster
C     appears or all the sites are filled.  This can be repeated for the
C     same pb value by adjusting the number of iterations.
C
C     Details for how the site-bond code works is given in the
C     appropriate bondsite.f file.
C
C     The output from the code will be the fraction of sites filled, the
C     fraction of bonds filled, and the random numbers used for that
C     particular iteration (for checking or debugging).  If the fraction
C     of sites filled is equal to zero, then that means that no lattice-
C     spanning cluster was present.
C
C     The critical percolation thresholds for an infinite square
C     lattice are 1/2 and 0.593 for the bond and site problems,
C     respectively.
C
C     Last Revision: November 19, 2010
C     ==================================================================

      program bs_perc

      implicit none
      integer m, n, t, cln, pbc, scn, bcn, rc, nb
      integer lcn, lcs, clsum, topfill, botfill, perccln, perccls
      integer maxcn, maxcs, i, j, k, l, temp, tbonds, rlc
      integer seed, iter, pbcount, ii, jj, kk
      integer s(1000000), sorder(1000000), c(1000000), nn(10)
      integer b(3000000,3), border(3000000,2), btemp(2), nnb(12,4)
      integer sseed(1000), bseed(1000), pseed(100)
      real rand
      double precision ps, pb, pbarray(100)
      external nearestn
      common m, n, t, pbc, nn, scn

      open(unit=10, file='bs_perc.txt')

C     dimensions of lattice m x n

         m = 10                   !width
         n = 10                   !height
         t = m*n                  !total number of lattice sites

C     main random number seed for generating random numbers in arrays
C     this seed is the only user-specified seed

         seed = 229102

C     toggle for periodic boundary conditions on the left/right sides

         pbc = 0                  !set to 1 if using p.b.c.

C     specify the array of pb values to use

         do i = 1,100             !zero out arrays
            pbarray(i) = 0
            pseed(i) = 0
         end do
         pbcount = 71             !number of pb 'points'
C        constants in the following loop can need be adjusted depending
C        upon the discretization desired
         do i = 1, pbcount
            pbarray(i) = 0.30d+00+(0.01d+00*(i-1))
         end do

C     specify the number of iterations per pb point

         iter = 1000

C     lattice-specific geometry constants (need to be changed if the
C     lattice is changed)

         scn = 4                  !site coordination number
         bcn = 6                  !bond coordination number

C     calculate the overall total number of bonds possible (depends on
C     lattice geometry

         if (pbc.eq.1) then
            nb = m*(2*n-1)
         else
            nb = (2*m*n)-m-n
         end if

C     populate pseed array with random numbers

         call srand(seed)
         do i=1,100               !second number must be equal to declared array size
            pseed(i) = int(rand(0)*10000000)+1
         end do


C     ******************************************************************
C     begin main code loop
C     ******************************************************************

      do ii = 1,pbcount

C     populate the site and bond random number arrays

         call srand(pseed(ii))    !seeds random number generator
         do jj = 1,1000           !second number must be equal to declared array size
            sseed(jj) = int(rand(0)*10000000)+1
            bseed(jj) = int(rand(0)*10000000)+1
         end do

C     run site-bond problems to percolation

         do jj = 1,iter

C     ------------------------------------------------------------------
C     begin (edited) site-bond code
C     ------------------------------------------------------------------

C     initialize arrays to zero for all elements (if the declared array
C     sizes change, you will need to modify the do loops here - the
C     bond array should have 3x the number of rows as the site array to
C     account for the geometry)

         do i = 1,1000000
            s(i) = 0              !site array
            c(i) = 0              !cluster array
         end do

         do i = 1,10
            nn(i) = 0             !nearest neighbor array
         end do

         do i = 1,3000000
            do j = 1,3
               b(i,j)=0           !bond array
            end do
         end do

C     initialize constants for checking for the infinite cluster

         maxcn = 0                !number of the largest overall cluster
         maxcs = 0                !size of the largest overall cluster
         perccln = 0              !lattice-spanning cluster number
         perccls = 0              !lattice-spanning cluster size

C     generate the array that lists the order of site occupation

         call srand(sseed(jj))    !seeds the random number generator

C        initialize sorder array

         do i = 1,t
            sorder(i) = i
         end do

C        permute sites randomly in the sorder array

         do i = 1,t
            j = i + (t-i+1)*rand(0)
            temp = sorder(i)
            sorder(i) = sorder(j)
            sorder(j) = temp
         end do

C     generate the list of bonds in the bond array and set the
C     occupation equal to zero; initialize the bond order array

         rc = 1                   !row number (counter)

         do i = 1,(t-1)
            call nearestn(i)
            do j = 1, scn
               if (nn(j).gt.i) then
                  b(rc,1) = i
                  b(rc,2) = nn(j)
                  b(rc,3) = 0
                  border(rc,1) = i
                  border(rc,2) = nn(j)
               else
                  goto 210
               end if
               rc = rc + 1
210            continue
            end do
         end do

C     randomize the order in which the bonds are filled in the bond
C     order array

         call srand(bseed(jj))    !seeds the random number generator

C     permute bonds randomly in the bond order array

         do i = 1,nb
            j = i + (nb-i+1)*rand(0)
            btemp(1) = border(i,1)
            btemp(2) = border(i,2)
            border(i,1) = border(j,1)
            border(i,2) = border(j,2)
            border(j,1) = btemp(1)
            border(j,2) = btemp(2)
         end do

C     occupy the specified fraction of bonds

         tbonds = pbarray(ii)*nb  !number of bonds to occupy
         cln = 1                  !lowest unused cluster number

         do i = 1,tbonds
            do j=1,nb             !this loop occupies the bond
               if (b(j,1).eq.border(i,1)) then
                  if (b(j,2).eq.border(i,2)) then
                     b(j,3) = cln
                  end if
               end if
            end do
            c(cln) = 1            !all bonds occupy clusters of size 1
            cln = cln+1           !increment the lowest unused cluster number
         end do

         pb = real(tbonds)/real(nb)

C     set the maximum cluster size to 1 and arbitrarily set the largest
C     cluster number to be 1

         maxcn = 1
         maxcs = 1

C     occupy sites and update cluster numbers/sizes appropriately

         lcn = 0                  !temp variable for largest cluster number
         lcs = 0                  !temp variable for largest cluster size

      do i = 1,t
C     initialize the nearest neighbor bond matrix to zero

            do j = 1,4
               do k = 1,12
                  nnb(k,j) = 0
               end do
            end do

C     populate the nearest neighbor bond matrix with useful info

            rc = 1                !reset row number/counter

            call nearestn(sorder(i))
            do j = 1,scn
               if (nn(j).ne.0) then
                  if(nn(j).gt.sorder(i)) then
                     nnb(rc,1) = sorder(i)
                     nnb(rc,2) = nn(j)
                  else
                     nnb(rc,1) = nn(j)
                     nnb(rc,2) = sorder(i)
                  end if
                  rc = rc + 1
               end if
            end do

            do j = 1,nb
               do k = 1,scn
                  if(b(j,1).eq.nnb(k,1)) then
                     if(b(j,2).eq.nnb(k,2)) then
                        nnb(k,3) = b(j,3)        !cluster number
                        nnb(k,4) = c(b(j,3))     !cluster size
                     end if
                  end if
               end do
            end do

C     find out which nearest neighbor bond is attached to the largest
C     cluster and get the largest cluster number and size associated
C     with that nearest neighbor bond

            lcn = nnb(1,3)        !current largest cluster number (arbitrarily set)
            lcs = nnb(1,4)        !current largest cluster size (arbitrarily set)
            rlc = 1               !row of data in the n.n. bond matrix containing the largest cluster

            do k = 2,scn
               if (nnb(k,1).ne.0) then           !no valid n.n. bond if =0
                  if (nnb(k,3).ne.0) then        !n.n. bond unoccupied if =0
                     if (nnb(k,4).gt.lcs) then
                        lcn = nnb(k,3)
                        lcs = nnb(k,4)           !or lcs = c(nnb(k,3))
                        rlc = k
                     end if
                  end if
               end if
            end do


C     Case 1: no n.n. bonds are occupied; site becomes its own cluster

            if (lcs.eq.0) then
               s(sorder(i)) = cln
               c(cln) = 1
               cln = cln+1
               goto 310
            end if

C     Case 2: one or more nearest neighbors occupied (i.e. lcs > 0)

            clsum = lcs                                    !temp variable for summing cluster size

            do k = 1,scn
               if (nnb(k,1).ne.0) then                     !no valid n.n. if =0
                  if (nnb(k,3).ne.0) then                  !n.n. unoccupied if =0
                     if (nnb(k,3).ne.lcn) then             !filters out n.n. already belonging to the largest cluster
                        if (k.ne.1) then                   !this should check to make sure the same cluster number is not "counted" more than once
                           do l = 1,(k-1)
                              if (nnb(l,3).eq.nnb(k,3)) then
                                 goto 320
                              end if
                           end do
                        end if
                        clsum = clsum + nnb(k,4)           !adds cluster size to largest cluster
                        do j = 1,nb
                           if (b(j,3).eq.nnb(k,3)) then    !this loop re-numbers all the bonds belonging to the smaller cluster with the largest cluster number
                              b(j,3) = lcn
                           end if
                        end do
                        do j = 1,t                         !this loop re-numbers all the sites belonging to the smaller cluster with the largest cluster number
                           if (s(j).eq.nnb(k,3)) then
                              s(j) = lcn
                           end if
                        end do
320                     c(nnb(k,3)) = 0                    !smaller cluster now contains zero sites
                     end if
                  end if
               end if
            end do

C     add occupied site to the largest cluster

            s(sorder(i)) = lcn
            clsum = clsum+1
            c(lcn) = clsum

C     update the largest overall cluster number and size variables
310         if (c(lcn).gt.maxcs) then
               maxcs = c(lcn)
               maxcn = lcn
            end if

C     calculate the fraction of sites occupied

         ps = real(i)/real(t)

C     check for percolation, i.e. presence of lattice spanning cluster
         if (i.ge.n) then         !must occupy a minimum of n sites before possible percolation
         do k = 1,(cln-1)
            botfill = 0           !used for checking if the bottom row contains a site
            topfill = 0           !used for checking if the top row contains a site
            if (c(k).ge.(2*n-1)) then   !only clusters with size 2n-1 could span the lattice
               do kk = 1,m
                  if (s(kk).eq.k) then
                     botfill = 1
                     goto 410
                  end if
               end do
               goto 430
410            do kk = (t-m+1),t
                  if (s(kk).eq.k) then
                     topfill = 1
                     goto 420
                  end if
               end do
420            if ((topfill+botfill).eq.2) then
                  perccln = k     !infinite cluster number
                  perccls = c(k)  !infinite cluster size
                  goto 440
               end if
430         end if
         end do
         end if
      end do

C     ------------------------------------------------------------------
C     end site-bond code
C     ------------------------------------------------------------------

440      if (perccln.eq.0) then
            write(6,111) sseed(jj), bseed(jj), 0.00d+00, pb
            write(10,111) sseed(jj), bseed(jj), 0.00d+00, pb
         else
            write(6,111) sseed(jj), bseed(jj), ps, pb
            write(10,111) sseed(jj), bseed(jj), ps, pb
         end if

         end do   !'end do' for one value of pb

      end do

C     ******************************************************************
C     end main code loop
C     ******************************************************************

111      format(i10,",",i10,",",f12.9,",",f12.9)

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
