C     ==================================================================
C     Program: sb_perc.f
C     Author: Isaiah Steinke
C
C     Runs multiple site-bond problems on a square lattice.  The
C     user must specify an array of values for ps at which to run the
C     site-bond problem.  The code will occupy the specified fraction of
C     sites and then occupy bonds until a lattice-spanning cluster
C     appears or all the bonds are filled.  This can be repeated for the
C     same ps value by adjusting the number of iterations.
C
C     Details for how the site-bond code works is given in the
C     appropriate sitebond.f file.
C
C     The output from the code will be the fraction of sites filled, the
C     fraction of bonds filled, and the random numbers used for that
C     particular iteration (for checking or debugging).  If the fraction
C     of bonds filled is equal to zero, then that means that no lattice-
C     spanning cluster was present.
C
C     The critical percolation thresholds for an infinite square
C     lattice are 1/2 and 0.593 for the bond and site problems,
C     respectively.
C
C     Last Revision: November 19, 2010
C     ==================================================================

      program sb_perc

      implicit none
      integer m, n, t, cln, pbc, scn, bcn, rc, nb
      integer lcn, lcs, oldcn, clsum, topfill, botfill, perccln, perccls
      integer maxcn, maxcs, i, j, k, temp, tsites
      integer seed, iter, pscount, ii, jj, kk
      integer s(1000000), sorder(1000000), c(1000000), nn(10)
      integer b(3000000,3), border(3000000,2), btemp(2)
      integer sseed(1000), bseed(1000), pseed(100)
      real rand
      double precision ps, pb, psarray(100)
      external nearestn
      common m, n, t, pbc, nn, scn

      open(unit=10, file='sb_perc.txt')

C     dimensions of lattice m x n

         m = 50                   !width
         n = 50                   !height
         t = m*n                  !total number of lattice sites

C     main random number seed for generating random numbers in arrays
C     this seed is the only user-specified seed

         seed = 8811064

C     toggle for periodic boundary conditions on the left/right sides

         pbc = 0                  !set to 1 if using p.b.c.

C     specify the array of ps values to use

         do i = 1,100             !zero out arrays
            psarray(i) = 0
            pseed(i) = 0
         end do
         pscount = 42             !number of ps 'points'
C        constants in the following loop can need be adjusted depending
C        upon the discretization desired
         do i = 1, pscount
            psarray(i) = 0.59d+00+(0.01d+00*(i-1))
         end do

C     specify the number of iterations per ps point

         iter = 100

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

      do ii = 1,pscount

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

         do i=1,1000000
            s(i) = 0              !site array
            c(i) = 0              !cluster array
         end do

         do i=1,10
            nn(i) = 0             !nearest neighbor array
         end do

         do i=1,3000000
            do j=1,3
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

C        initialize order array

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

C     occupy the specified fraction of sites

         tsites = psarray(ii)*t   !number of sites to occupy
         cln = 1                  !lowest unused cluster number

         do i=1,tsites
            s(sorder(i)) = cln    !occupies site
            c(cln) = 1            !all sites occupied are of size 1
            cln = cln+1           !increment the lowest unused cluster number
         end do

         ps = real(tsites)/real(t)

C     set the maximum cluster size to 1 and arbitrarily set the largest
C     cluster number to be 1

         maxcn = 1
         maxcs = 1

C     occupy bonds and update cluster numbers/sizes appropriately

         lcn = 0                  !temp variable for largest cluster number
         lcs = 0                  !temp variable for largest cluster size

      do i=1,nb
C     check the two sites that the bond connects
         do j=1,nb
            if (b(j,1).eq.border(i,1)) then
               if (b(j,2).eq.border(i,2)) then
C     Case 1: neither site is occupied, bond becomes its own cluster
                  if (s(b(j,1)).eq.0) then
                     if (s(b(j,2)).eq.0) then
                        b(j,3) = cln
                        c(cln) = 1
                        cln = cln+1
                        goto 310
                     end if
                  end if
C     Case 2: only one site is occupied, bond is added to cluster of
C     occupied site
                  if (s(b(j,1)).gt.0) then
                     if (s(b(j,2)).eq.0) then
                        lcn = s(b(j,1))
                        lcs = c(s(b(j,1)))
                        b(j,3) = lcn
                        c(s(b(j,1))) = lcs+1
                        goto 310
                     end if
                  end if
                  if (s(b(j,1)).eq.0) then
                     if (s(b(j,2)).gt.0) then
                        lcn = s(b(j,2))
                        lcs = c(s(b(j,2)))
                        b(j,3) = lcn
                        c(s(b(j,2))) = lcs+1
                        goto 310
                     end if
                  end if
C     Case 3: both sites are occupied; add everything to the larger
C     cluster and zero out the smaller cluster
C        this loop checks to see if both sites are already part of the
C        same cluster
                  if (s(b(j,1)).eq.s(b(j,2))) then
                     lcn = s(b(j,1))
                     lcs = c(s(b(j,1)))
                     b(j,3) = lcn
                     c(s(b(j,1))) = lcs+1
                     goto 310
                  end if
C        this loop is for sites belonging to different clusters
                  if (c(s(b(j,1))).gt.c(s(b(j,2)))) then
                     lcn = s(b(j,1))
                     oldcn = s(b(j,2))
                     lcs = c(s(b(j,1)))
                     b(j,3) = lcn
                     clsum = lcs+c(oldcn)+1
                     do k=1,t
                        if (s(k).eq.oldcn) then
                           s(k) = lcn
                        end if
                     end do
                     do k=1,nb
                        if (b(k,3).eq.oldcn) then
                           b(k,3) = lcn
                        end if
                     end do
                     c(oldcn) = 0
                     c(s(b(j,1))) = clsum
                     goto 310
                  else
                     lcn = s(b(j,2))
                     oldcn = s(b(j,1))
                     lcs = c(s(b(j,2)))
                     b(j,3) = lcn
                     clsum = lcs+c(oldcn)+1
                     do k=1,t
                        if (s(k).eq.oldcn) then
                           s(k) = lcn
                        end if
                     end do
                     do k=1,nb
                        if (b(k,3).eq.oldcn) then
                           b(k,3) = lcn
                        end if
                     end do
                     c(oldcn) = 0
                     c(s(b(j,2))) = clsum
                     goto 310
                  end if
               end if
            end if
         end do
C     update the largest overall cluster number and size variables
310      if (c(lcn).gt.maxcs) then
            maxcs = c(lcn)
            maxcn = lcn
         end if

C     calculate the fraction of bonds occupied

         pb = real(i)/real(nb)

C     check for percolation, i.e. presence of lattice spanning cluster
         if (i.ge.(n-1)) then     !must occupy a minimum of n-1 bond before possible percolation
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
            write(6,111) sseed(jj), bseed(jj), ps, 0.00d+00
            write(10,111) sseed(jj), bseed(jj), ps, 0.00d+00
         else
            write(6,111) sseed(jj), bseed(jj), ps, pb
            write(10,111) sseed(jj), bseed(jj), ps, pb
         end if

         end do   !'end do' for one value of ps

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
