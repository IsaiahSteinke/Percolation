C     ==================================================================
C     Program: site_perc.f
C     Author: Isaiah Steinke
C
C     Runs multiple trials of the site percolation code on a
C     square lattice in two dimensions.  Details for the operation
C     of the basic site percolation code are given in the program
C     "site.f".  In this code, the user specifies the lattice geometry
C     and the number of trials to run.  Each trial is run until a
C     lattice-spanning (or "infinite cluster") is present.
C
C     The code outputs a single text file of data containing the
C     following:
C        -the random number used to seed each trial (for debugging
C         purposes).
C        -the fraction of sites filled when each trial was terminated
C        -the number of sites (i.e. the size) of the infinite cluster
C        -the size of the infinite cluster
C
C     *Critical percolation threshold for the site problem on an
C     infinite square lattice is 0.593.
C
C     Last Revision: July 8, 2010
C     ==================================================================

      program site_perc

      implicit none
      integer m, n, t, cln, pbc, sn, sf, scn, bcn, rc, nb
      integer seed, lcn, lcs, nnlc, oldcn, clsum, numtrials
      integer maxcn, maxcs, topfill, botfill, perccln, perccls
      integer i, j, k, l, ii, temp
      integer s(1000000), c(1000000), nn(10), order(1000000)
      integer tseed(50000)
      real rand
      double precision f
      external nearestn
      common m, n, t, pbc, nn, scn

      open(unit=10, file='site_perc.txt')

C     dimensions of lattice m x n

         m = 50                   !width
         n = 50                   !height
         t = m*n                  !total number of lattice sites

C     toggle for periodic boundary conditions on the left/right sides

         pbc = 0                  !set to 1 if using p.b.c.

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

C     number of percolation trials to execute and random number seed

         numtrials = 1000         !upper limit: declared array size for tseed
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

C     initialize arrays to zero for all elements

         do i=1,1000000        !second number should be equal to the declared array size
            s(i) = 0           !site array
            c(i) = 0           !cluster array
         end do

C     initialize constants for checking for the infinite cluster

         maxcn = 0                !number of the largest overall cluster
         maxcs = 0                !size of the largest overall cluster

C     ------------------------------------------------------------------
C     generate the array that lists the order of site occupation
C     ------------------------------------------------------------------

         call srand(tseed(ii))

C        initialize order array

         do i = 1,t
            order(i) = i
         end do

C        permute sites randomly in the order array

         do i = 1,t
            j = i + (t-i+1)*rand(0)
            temp = order(i)
            order(i) = order(j)
            order(j) = temp
         end do
C     ------------------------------------------------------------------

C     ------------------------------------------------------------------
C     fill lattice with sites to percolation
C     ------------------------------------------------------------------

C     initialize constants/counters/calculated values

         sf = 0                   !number of sites filled
         f = real(sf)/real(t)     !fraction of sites filled
         cln = 1                  !lowest unused cluster number

         do i = 1,t
            sn = order(i)

C     initialize constants for tracking largest cluster characteristics

            lcn = 0               !cluster number of largest cluster
            lcs = 0               !size of largest cluster
            nnlc = 0              !n.n. with largest cluster

            call nearestn(sn)

C     find the nearest neighbor with the largest cluster size and get
C     the largest cluster number and cluster size associated with that
C     nearest neighbor

            lcn = s(nn(1))
            lcs = c(s(nn(1)))
            nnlc = nn(1)

            do k = 2,scn
               if (nn(k).ne.0) then              !no valid n.n. if =0
                  if(s(nn(k)).ne.0) then         !n.n. unoccupied if =0
                     if (c(s(nn(k))).gt.lcs) then
                        lcn = s(nn(k))
                        lcs = c(s(nn(k)))
                        nnlc = nn(k)
                     end if
                  end if
               end if
            end do

C     Case 1: no nearest neighbors occupied (i.e. lcs = 0)

            if (lcs.eq.0) then
               s(sn) = cln
               c(cln) = 1
               cln = cln+1
               goto 220
            end if

C     Case 2: one or more nearest neighbors occupied (i.e. lcs > 0)

            clsum = lcs                               !temp variable for summing cluster size

            do k = 1,scn
               if (nn(k).ne.0) then                   !no valid n.n. if =0 (i.e. corner/side site)
                  if (s(nn(k)).ne.0) then             !n.n. unoccupied if =0
                     if (s(nn(k)).ne.s(nnlc)) then    !filters out n.n. already belonging to the largest cluster
                        if (k.ne.1) then              !this should check to make sure the same cluster number
                           do l = 1, (k-1)            !is not "counted" more than once
                              if(s(nn(l)).eq.s(nn(k))) then
                                 goto 260
                              end if
                           end do
                        end if
                        clsum = clsum + c(s(nn(k)))   !adds cluster size to largest cluster
                        oldcn = s(nn(k))              !temp variable for the old cluster number
                        do j = 1,t
                           if (s(j).eq.oldcn) then    !this loop re-numbers all the sites belonging to the smaller
                              s(j) = lcn              !cluster with the largest cluster number
                           end if
                        end do
260                     c(oldcn) = 0                  !smaller cluster now contains zero sites
                     end if
                  end if
               end if
            end do

C     add randomly picked site to the largest cluster

            s(sn) = lcn
            clsum = clsum+1
            c(lcn) = clsum

220         sf = sf + 1
            f = real(sf)/real(t)

C     this small loop should check for the largest overall cluster
C     and update the variables for the overall largest cluster size
C     and number

            if (c(lcn).gt.maxcs) then
               maxcs = c(lcn)
               maxcn = lcn
C           this sub-loop fixes the case where only one site is occupied
            else if (lcs.eq.0) then
               if (maxcs.eq.0) then
                  maxcs = 1
                  maxcn = 1
               end if
            end if

C     check for percolation, i.e. presence of lattice spanning cluster

         if (i.ge.n) then         !need to occupy at least n sites before infinite cluster is possible
            do k = 1,(cln-1)
               botfill = 0        !used for checking if the bottom row contains a site
               topfill = 0        !used for checking if the top row contains a site
               if (c(k).ge.n) then   !only clusters with size of at least the height of the lattice could span the lattice
                  do j = 1,m
                     if (s(j).eq.k) then
                        botfill = 1
                        goto 230
                     end if
                  end do
                  goto 270
230               do j = (t-m+1),t
                     if (s(j).eq.k) then
                        topfill = 1
                        goto 240
                     end if
                  end do
240               if ((topfill+botfill).eq.2) then
                     perccln = k      !infinite cluster number
                     perccls = c(k)   !infinite cluster size
                     goto 250
                  end if
270            end if
            end do
         end if

         end do

C     output information when infinite cluster is present
250      write(6,*) tseed(ii), f, maxcs, perccls
         write(10,111) tseed(ii), f, maxcs, perccls

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
