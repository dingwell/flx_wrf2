      subroutine calcprec(pacc,prate)
C                           i    o
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

      real pacc(0:nxmax-1,0:nymax-1,1,2)
      real prate(0:nxmax-1,0:nymax-1,1,2)
      real dt,smallnum

      integer ix,jy

C     data smallnum /1.e-38/ ! Some indata had values of ~ -10⁻⁷
      data smallnum /1.e-07/

C Get dt in hours:
      dt = real(memtime(2)-memtime(1))/3600.
      !write(*,*) 'AD check dt:',dt
      !write(*,*) 'AD check memtime:',memtime(1),memtime(2)

C Get precipitation rate 
      do jy=0,nymin1
      do ix=0,nxmin1
        ! Check for 3 possible cases:
        ! 1: prate .ge. 0:
        if ( pacc(ix,jy,1,memind(2))
     +  .ge. pacc(ix,jy,1,memind(1)) ) then
          prate(ix,jy,1,memind(2)) =
     +    ( pacc(ix,jy,1,memind(2))-pacc(ix,jy,1,memind(1)) )/dt

        ! 2: prate .lt. 0 due to numeric error -> set to zero:
        else if ( pacc(ix,jy,1,memind(2))
     +  .gt. pacc(ix,jy,1,memind(1))-smallnum*dt) then
          prate(ix,jy,1,memind(2)) = 0

        ! 3: prate .lt. 0, probably due to error in input files
        else
          write(*,9100) 'ERROR: Decrease in accumulated precipitation.'
          write(*,9110) 'Possible cold-start of WRF-model encountered,'
          write(*,9110) 'unable to calculate precipitation rate'
          write(*,9120) pacc(ix,jy,1,memind(1)),pacc(ix,jy,1,memind(2))
          write(*,9130) memtime(memind(1)),memtime(memind(2))
          write(*,9140) ix,jy
          stop
        endif

      end do
      end do

9100  format( / '*** calcprec          -- ', a )
9110  format(   '***                      ', a )
9120  format(   '*** pacc(prev,next) =    ', 2(1x,e14.6) )
9130  format(   '*** memtime(prev,next) = ', 2(1x,i8.1) )
9140  format(   '*** ix, jy =             ', 2(1x,i3) )

      end
