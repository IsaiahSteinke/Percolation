C     ==================================================================
C     Program: permute.f
C     Author: Isaiah Steinke
C
C     This program takes an array of ordered integer numbers starting
C     from 1 up to t (the total specified number of elements) that are
C     incremented by 1 and randomly permutes the order of integers.
C
C     The idea for this is adapted from the C code presented in:
C     Newman & Ziff, Phys. Rev. E, Vol. 64, 016706 (2001).
C
C     Last Revision: November 5, 2009
C     ==================================================================

      implicit none
      integer i, j, temp, t, seed
      integer order(1000000)
      real rand
      
      open(unit=10, name='permute.txt')

C     specify seed for random number generator and seed

      seed = 4711904
      call srand(seed)

C     total number of elements in the array

      t = 10000
      
C     initialize order array

      do i = 1, t
         order(i) = i
      end do
      
C     permute sites randomly in the order array

      do i = 1, t
         j = i + (t-i+1)*rand(0)
         temp = order(i)
         order(i) = order(j)
         order(j) = temp
      end do
      
C     output array of permuted integers

      do i = 1, t
         write(6,*) order(i)
         write(10,*) order(i)
      end do
      
      end
