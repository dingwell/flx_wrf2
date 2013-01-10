      subroutine calcprec_nests(paccn,praten)
C                                  i    o
***********************************************************************
*                                                                     *
* This subroutine calculates hourly precipitation rate based on       *
* accumulated precipitation from WRF.                                 *
*                                                                     *
***********************************************************************
*                                                                     *
* AUTHOR:      A. DINGWELL                                            *
* DATE:        2012-12-13                                             *
*                                                                     *
***********************************************************************
*                                                                     *
* VARIABLES:                                                          *
* pacc  [mm]          Accumulated precipitation                       *
* prate [mm/h]        precipitation rate                              *
*                                                                     *
* dt [h]              time between the WRF outputs                    *
*                                                                     *
* smallnum            Lowest handled precipitation rate [mm/h]        *
*                                                                     *
***********************************************************************

      include 'includepar'
      include 'includecom'

      real paccn(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
      real praten(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
      real dt,smallnum

      integer ix,jy,l

C     data smallnum /1.e-38/ ! Some indata had values of ~ -10⁻⁷
      data smallnum /1.e-07/

C Get dt in hours:
      dt = real(memtime(2)-memtime(1))/3600.

C Check that acc precipitation increases with time and calc rate
      do l=1,numbnests
      do jy=0,nyn(l)-1
      do ix=0,nxn(l)-1
        ! Check for 3 possible cases:
        ! 1: prate .ge. 0
        if ( paccn(ix,jy,1,memind(2),l)
     +  .ge. paccn(ix,jy,1,memind(1),l) ) then
          praten(ix,jy,1,memind(2),l) =
     +    (paccn(ix,jy,1,memind(2),l)-paccn(ix,jy,1,memind(1),l) )/dt
          
        ! 2: prate .lt. 0 due to numeric error -> set to zero
        else if ( paccn(ix,jy,1,memind(2),l)
     +  .gt. paccn(ix,jy,1,memind(1),l)-smallnum*dt) then
          praten(ix,jy,1,memind(2),l) = 0

        ! 3: prate .lt. 0, probably due to error in input files
        else
          write(*,9100) 'ERROR: Decrease in accumulated precipitation.'
          write(*,9110) 'Possible cold-start of WRF-model encountered,'
          write(*,9110) 'unable to calculate precipitation rate'
          write(*,9120) paccn(ix,jy,1,memind(1),l),
     +      paccn(ix,jy,1,memind(2),l)
          write(*,9130) memtime(memind(1)),memtime(memind(2))
          write(*,9140) ix,jy,l
          stop
        endif
      end do
      end do
      end do

9100  format( / '*** calcprec_nests    -- ', a )
9110  format(   '***                      ', a )
9120  format(   '*** pacc(prev,next) =    ', 2(1x,e14.6) )
9130  format(   '*** memtime(prev,next) = ', 2(1x,i8.1) )
9140  format(   '*** ix, jy =             ', 2(1x,i3) )

      end
