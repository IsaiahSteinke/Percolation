C     ==================================================================
C     Program: randtest.f
C     Author: Isaiah Steinke
C
C     A short program to test and illustrate the use of the random
C     number generator in Fortran, rand.
C
C     *The function rand is limited to only real number output.*
C     It appears that this limits the output to 9 significant figures.
C
C     Last Updated: September 25, 2009
C     ==================================================================

      implicit none
      integer seed, i
      real rand
      integer time(3)

C     seed the random number generator (can be user specified).
         seed = 73789983               !must be an integer of kind 4
         call srand(seed)

C     print the first 10 numbers of the sequence.
C     rand(0) returns the next number in the sequence.

         do i = 1, 10
            write(6,*) rand(0)
         end do

         write(6,*)

C     re-seed and write out the first 10 numbers as a check to see if
C     the sequence can be restarted (output should be identical to the
C     first sequence of 10 numbers).

         call srand(seed)
         do i = 1, 10
            write(6,*) rand(0)
         end do
         
         write(6,*)
         
C     re-seed and write out the first 10 numbers as integers between
C     1 and 100.

         call srand(seed)
         do i = 1, 10
            write(6,*) int(rand(0)*100)+1
         end do
         
         write(6,*)
         
C     re-seed with current time and write out first 10 numbers.
C     these should change with execution of the program.

         call itime(time)              !gets the system time
         seed = time(1)+time(2)+time(3)
         call srand(seed)
         
         do i = 1, 10
            write(6,*) rand(0)
         end do

         write(6,*)
         
C     re-seed via the rand function and write out the first 10 numbers.
C     as long as the argument for rand is not either 0 or 1, the integer
C     should be used as the new seed.

         write(6,*) rand(34543234)     !must be an integer
         do i = 1, 9
            write(6,*) rand(0)
         end do

      end
