C     ==================================================================
C     Program: bond_cond.f
C     Author: Isaiah Steinke
C
C     Runs a bond percolation problem on a two dimensional triangular
C     lattice and calculates the conductance of the lattice-spanning
C     cluster as bonds are added. This program is set up to run multiple
C     instances of the bond percolation problem.  Details of the bond
C     percolation code are given in either the bond.f or bondc.f files.
C
C     The output from this code primarily consists of the fraction of
C     bonds occupied (pb), the calculated conductances (Gbot and Gtop),
C     and the average of the two conductances for each realization
C     of the bond problem run.
C
C     *Critical percolation threshold for the bond problem on an
C     infinite triangular lattice is 0.3473.
C
C     This code uses subroutines to calculate the conductance from:
C     "Numerical Recipes in Fortran 77: The Art of Scientific
C     Computing", by William Press, et. al., Cambridge Press, 2001.
C     http://www.nr.com
C
C     Last Revision: October 28, 2010
C     ==================================================================

      program bond_cond

      implicit none
      integer m, n, t, pbc, rc, nb, bf, topfill, botfill
      integer scn, bcn, cln, lcn, lcs, rlc, clsum, numtrials
      integer i, j, k, l, seed, perccln, ii, jj, tseed(1000)
      integer b(1000000,3), border(1000000,2), c(1000000), nbarr(250)
      integer nn(10), btemp(2), nnb(25,4), ija(20000), iter
      real rand
      double precision pb, pc, pbarr(250)
      double precision Va, g0, G(2500,2500), Itemp(2500), Iout(2500)
      double precision rowsum, Gtemp(2500,2500), Vint(2500)
      double precision V(2500), Itop, Ibot, Gtop, Gbot
      double precision sa(20000), err
      external nearestn, snrm
      common m, n, t, pbc, nn, scn
      COMMON /mat/ sa, ija

      open(unit=10, file='bondcond.txt')

C     lattice dimensions, m x n

         m = 10                   !number of lattice sites wide
         n = 10                   !number of lattice sites high
         t = m*n                  !total no. of lattice sites

C     toggle for periodic boundary conditions on the left/right sides
C     note that the pbc flag only makes sense for m even with the
C     current lattice geometry (triangular)

         pbc = 0                  !set to 1 if using p.b.c.

C     Conductance calculation parameters

         Va = 1.00d+00            !applied voltage [V]
         g0 = 1.00d+00            !base bond conductance [1/ohm]
         
C     number of percolation trials to execute and random number seed

         numtrials = 1            !upper limit: declared array size for tseed
         seed = 58302             !must be an integer
         call srand(seed)         !seed random number generator
C        generate the array of random number seeds for each trial
         do i=1,1000              !second number must be equal to declared array size
            tseed(i) = int(rand(0)*10000000)+1
         end do

C     calculate the overall total number of bonds possible (depends on
C     lattice geometry)

         if (pbc.eq.1) then
            nb = m*(3*n-2)
         else
            nb = (3*m*n)-(2*m)-(2*n)+1
         end if
         
C     array of pb values at which to calculate the conductance (may
C     need to be adjusted based on lattice size)

         do i = 1, 250
            pbarr(i) = 0.00d+00
            nbarr(i) = 0
         end do

         pbarr(1) = 0.35d+00
         do i = 2,131
            pbarr(i) = pbarr(i-1)+5.00d-03
         end do
         do i = 1,250
            nbarr(i) = pbarr(i)*nb
C            write(6,*) pbarr(i), nbarr(i)
C            write(10,*) pbarr(i), nbarr(i)
         end do

C     lattice-specific geometry constants (need to be changed if the
C     lattice is changed)

         scn = 6                  !site coordination number
         bcn = 10                 !bond coordination number

C     output setup parameters to file

         write(10,*) "m =", m
         write(10,*) "n =", n
         write(10,*) "t =", t
         write(10,*) "pbc =", pbc
         write(10,*) "total number of bonds in the lattice:", nb
         write(10,*) "Va =", Va
         write(10,*) "g0 =", g0
         write(10,*) "total number of iterations:", numtrials
         write(10,*) "random number for generating the trial seeds:"
     &   , seed
         write(10,*) "------------------------------"

C     ******************************************************************
C     begin main code loop
C     ******************************************************************

      do ii = 1,numtrials
         write(6,*) "Trial #", ii
         write(10,*) "Trial #", ii
         write(6,*) "Random number seed:", tseed(ii)
         write(10,*) "Random number seed:", tseed(ii)
         
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

C     ==================================================================
C     begin percolation code
C     ==================================================================

C     initialize constants/counters/calculated values

         bf = 0                   !number of bonds filled
         cln = 1                  !lowest unused cluster number
         jj = 1                   !counter for nbarr
         perccln = 0              !lattice-spanning cluster number
         
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
            pb = real(bf)/real(nb)

C     check for percolation, i.e. presence of lattice spanning cluster
         do l = 1,(cln-1)
            botfill = 0           !used for checking if the bottom row is connected to the cluster
            topfill = 0           !used for checking if the top row is connected to the cluster
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
340            do j = (t-m+1),t
                  do k = 1,nb
                     if (b(k,2).eq.j) then
                        if (b(k,3).eq.l) then
                           topfill = 1
                           goto 350
                        end if
                     end if
                  end do
               end do
350            if ((topfill+botfill).eq.2) then
                  if (perccln.eq.0) then
                     perccln = l     !infinite cluster number
                     pc = real(bf)/real(nb)
                     goto 370
                  else  !this part updates perccln in case it changes (gets absorbed into the max cluster)
                     perccln = l
                     goto 370
                  end if
               end if
360         end if
         end do

C     calculate the conductance of the infinite cluster (if present)
370      if (bf.eq.nbarr(jj)) then
         if (perccln.gt.0) then
            do l = 1,2500       !zero out conductance matrices
               do j = 1,2500
                  G(l,j) = 0.00d+00
                  Gtemp(l,j) = 0.00d+00
               end do
            end do
            do l = 1,2500       !zero out vectors used during calcuation
               Itemp(l) = 0.00d+00
               Vint(l) = 0.00d+00
               Iout(l) = 0.00d+00
               V(l) = 0.00d+00
            end do
            Itop = 0.00d+00
            Ibot = 0.00d+00
            Gtop = 0.00d+00
            Gbot = 0.00d+00
            do l = 1,nb           !populate off-diagonal elements of cond. matrix
               if (b(l,3).eq.perccln) then
                  G(b(l,1),b(l,2)) = -g0
                  G(b(l,2),b(l,1)) = G(b(l,1),b(l,2))
               else
                  G(b(l,1),b(l,2)) = -1.00d-12
                  G(b(l,2),b(l,1)) = G(b(l,1),b(l,2))
               end if
               if (b(l,1).gt.(t-2*m)) then !populate "temporary current" vector
                  if (b(l,1).le.(t-m)) then
                     if (b(l,2).gt.(t-m)) then
                        Itemp((b(l,1)-m)) = Itemp((b(l,1)-m))
     &                  -(G(b(l,1),b(l,2))*Va)
                     end if
                  end if
               end if
            end do
            do l = 1,t            !populate diagonal elements of cond. matrix
               rowsum = 0
               do j = 1,t
                  rowsum = rowsum + G(l,j)
               end do
               G(l,l) = -rowsum
            end do
            do l = 1,(t-2*m)      !populate Gtemp matrix properly
               do j = 1,(t-2*m)
                  Gtemp(l,j) = G((l+m),(j+m))
               end do
            end do
C           convert Gtemp matrix to row-indexed sparse storage mode
            do l = 1,20000
               sa(l) = 0.00d+00
               ija(l) = 0
            end do
            call sprsin(Gtemp,(t-2*m),2500,1.00d-16,20000,sa,ija)
C           set constants for linbcg routine
            iter = 0
            err = 1.00d-02
C           calculate internal node voltages using linbcg
            call linbcg((t-2*m),Itemp,Vint,2,1.00d-08,2500,iter,err)
            do l = 1,t            !populate voltage vector with the correct voltages
               if (l.le.m) then
                  V(l) = 0.00d+00
               else
                  if (l.gt.(t-m)) then
                     V(l) = Va
                  else
                     V(l) = Vint((l-m))  !Vint contains the solutions from linbcg routine
                  end if
               end if
            end do
C           convert G matrix to row-indexed sparse storage mode
            do l = 1,20000
               sa(l) = 0.00d+00
               ija(l) = 0
            end do
            call sprsin(G,t,2500,1.00d-10,20000,sa,ija)
C           calculate currents by multiplying G*V
            call dsprsax(sa,ija,V,Iout,t)
            do l = 1,m
               Ibot = Ibot + Iout(l)  !Iout contains solution from dsprsax routine
               Itop = Itop + Iout((l+t-m))
            end do
            Gtop = Itop/Va
            Gbot = abs(Ibot)/Va
            jj = jj+1
         else
            Gbot = 0.00d+00
            Gtop = 0.00d+00
            jj = jj+1
         end if
            write(6,111) pb, Gbot, Gtop, ((Gbot+Gtop)/2)
            write(10,111) pb, Gbot, Gtop, ((Gbot+Gtop)/2)
         end if
         
         end do

C     ==================================================================
C     end percolation code
C     ==================================================================

         write(6,*) "lattice-spanning cluster:", perccln
         write(10,*) "lattice-spanning cluster:", perccln
         write(6,*) "pc =", pc
         write(10,*) "pc =", pc
         write(6,*) "------------------------------"
         write(10,*) "------------------------------"

      end do

C     ******************************************************************
C     end main code loop
C     ******************************************************************
      

111   format(f12.9,",",f12.9,",",f12.9,",",f12.9)

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
          snrm=snrm+(sx(i)*sx(i))
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

