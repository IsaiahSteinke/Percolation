C     ==================================================================
C     Program: site.f
C     Author: Isaiah Steinke
C
C     Runs a site percolation problem on a square lattice in two
C     dimensions.  Lattice sites are randomly filled from the empty
C     lattice state up to a specified fraction of the lattice, ps.  It
C     is assumed that bonds are present between two nearest neighbor
C     lattice sites that are occupied.  Code utilizes two arrays for
C     percolation: the site array which lists which cluster number the
C     site belongs to, and the cluster array which lists the size of
C     the cluster.  Another array is used to determine the order the
C     lattice sites are filled by randomly swapping the ordered sites.
C
C     Code outputs four text files:
C        -siteorder.txt lists the order in which the lattice sites are
C         filled
C        -site.txt lists the site (or cluster) number, site array, and
C         cluster array
C        -siteocc.txt lists info for what happened when each site was
C         occupied (useful for debugging/checking)
C        -bondlist.txt lists all the bonds present in the lattice.
C
C     *Critical percolation threshold for the site problem on an
C     infinite square lattice is 0.593.
C
C     Last Revision: July 8, 2010
C     ==================================================================

      program site
      
      implicit none
      integer m, n, t, cln, pbc, sn, sf, scn, bcn, rc, nb
      integer seed, lcn, lcs, nnlc, oldcn, clsum
      integer maxcn, maxcs, topfill, botfill, perccln, perccls
      integer i, j, k, l, temp, tsites
      integer s(1000000), c(1000000), nn(10), order(1000000)
      integer blist(1000000,2)
      real rand
      double precision ps, f
      external nearestn
      common m, n, t, pbc, nn, scn
      
      open(unit=10, file='site.txt')
      open(unit=11, file='siteocc.txt')
      open(unit=12, file='siteorder.txt')
      open(unit=13, file='bondlist.txt')
      
C     dimensions of lattice m x n

         m = 50                   !width
         n = 50                   !height
         t = m*n                  !total number of lattice sites

C     toggle for periodic boundary conditions on the left/right sides

         pbc = 0                  !set to 1 if using p.b.c.

C     fraction of lattice sites filled (user specified)

         ps = 0.60d+00

C     lattice-specific geometry constants (need to be changed if the
C     lattice is changed)

         scn = 4                  !site coordination number
         bcn = 6                  !bond coordination number
         
C     initialize arrays to zero for all elements

         do i=1,1000000           !second number should be equal to the declared array size
            s(i) = 0              !site array
            c(i) = 0              !cluster array
         end do

         do i=1,10
            nn(i) = 0             !nearest neighbor array
         end do

         do i=1,1000000
            do j=1,2
               blist(i,j) = 0     !lists all the bonds in the lattice
            end do
         end do

C     calculate the overall total number of bonds possible (depends on
C     lattice geometry

         if (pbc.eq.1) then
            nb = m*(2*n-1)
         else
            nb = (2*m*n)-m-n
         end if

C     initialize constants for checking for the infinite cluster

         maxcn = 0                !number of the largest overall cluster
         maxcs = 0                !size of the largest overall cluster
         
C     ------------------------------------------------------------------
C     generate a list of all the possible bonds in the lattice (bonds
C     are listed by the sites they connect), output to text file,
C     bondlist.txt, for use in Matlab plotting routine
C     ------------------------------------------------------------------

         rc = 1                   !row number (counter)

         do i = 1, (t-1)
            call nearestn(i)
            do j = 1, scn
               if (nn(j).gt.i) then
                  blist(rc,1) = i
                  blist(rc,2) = nn(j)
               else
                  goto 310
               end if
               rc = rc + 1
310            continue
            end do
         end do
         
         do i = 1, nb
            write(13,121) blist(i,1), blist(i,2)
         end do
C     ------------------------------------------------------------------
         
C     ------------------------------------------------------------------
C     generate the array that lists the order of site occupation
C     ------------------------------------------------------------------

         seed = 1080115           !user defined seed, must be an integer
         call srand(seed)         !seeds the random number generator
         
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

C        output array of permuted integers

         do i = 1,t
            write(12,*) order(i)
         end do
C     ------------------------------------------------------------------

C     ==================================================================
C     main program loop
C     ==================================================================

C     initialize constants/counters/calculated values

         sf = 0                   !number of sites filled
         f = real(sf)/real(t)     !fraction of sites filled
         tsites = ps*t            !total number of sites to fill
         cln = 1                  !lowest unused cluster number

         do i = 1,tsites
            sn = order(i)
            write(6,*) "site chosen:", sn
            write(11,*) "site chosen:", sn

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
            write(6,*) "nearest neighbors:", (nn(j), j = 1,scn)
            write(11,*) "nearest neighbors:", (nn(j), j = 1,scn)

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
            
            write(6,*) "n.n. with largest cluster:", nnlc
            write(6,*) "largest cluster number:", lcn
            write(6,*) "largest neighbor cluster size:", lcs
            write(11,*) "n.n. with largest cluster:", nnlc
            write(11,*) "largest cluster number:", lcn
            write(11,*) "largest neighbor cluster size:", lcs

C     Case 1: no nearest neighbors occupied (i.e. lcs = 0)

            if (lcs.eq.0) then
               s(sn) = cln
               c(cln) = 1
               write(6,*) "*no n.n. occupied*"
               write(6,*) "site assigned to cluster number", cln
               write(11,*) "*no n.n. occupied*"
               write(11,*) "site assigned to cluster number", cln
               cln = cln+1
               goto 220
            end if

C     Case 2: one or more nearest neighbors occupied (i.e. lcs > 0)

            clsum = lcs                               !temp variable for summing cluster size
            write(6,*) "*one or more n.n. occupied*"
            write(11,*) "*one or more n.n. occupied*"

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
                write(6,*) "adding", c(s(nn(k))), " to largest cluster"
                write(6,*) "largest cluster is now", clsum
                write(11,*) "adding", c(s(nn(k))), " to largest cluster"
                write(11,*) "largest cluster is now", clsum
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
            write(6,*) "site assigned to cluster number", lcn
            write(6,*) "size of cluster number", lcn, " is now", c(lcn)
            write(11,*) "site assigned to cluster number", lcn
            write(11,*) "size of cluster number", lcn, " is now", c(lcn)
            
220         sf = sf + 1
            f = real(sf)/real(t)

            write(6,*) "fraction of lattice filled:", f
            write(6,*) "--------------------"
            write(11,*) "fraction of lattice filled:", f
            write(11,*) "--------------------"

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
            botfill = 0           !used for checking if the bottom row contains a site
            topfill = 0           !used for checking if the top row contains a site
            if (c(i).ge.n) then   !only clusters with size of at least the height of the lattice could span the lattice
               write(6,*) "testing cluster"   !useful debugging message
               write(11,*) "testing cluster"
               do j = 1,m
                  if (s(j).eq.i) then
                     botfill = 1
                     goto 230
                  end if
               end do
               write(6,*) "source end not connected"
               write(11,*) "source end not connected"
               goto 270
230            do j = (t-m+1),t
                  if (s(j).eq.i) then
                     topfill = 1
                     goto 240
                  end if
               end do
               write(6,*) "drain end not connected"
               write(11,*) "drain end not connected"
240            if ((topfill+botfill).eq.2) then
                  write(6,*) "infinite cluster present"
                  write(11,*) "infinite cluster present"
                  perccln = i     !infinite cluster number
                  perccls = c(i)  !infinite cluster size
                  write(6,*) "infinite cluster number:", perccln
                  write(11,*) "infinite cluster number:", perccln
                  write(6,*) "infinite cluster size:", perccls
                  write(11,*) "infinite cluster size:", perccls
                  goto 250
               end if
270         end if
         end do

         write(6,*) "no infinite cluster present"
         write(11,*) "no infinite cluster present"

250      write(6,*) "******************************"
         write(11,*) "******************************"

C     output site and cluster arrays

         do j = 1,t
            write(10,111) j, s(j), c(j)
         end do

111      format(i10,",",i10,",",i10)
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
