      subroutine check_ncerror(errcode)
C                             (integer)
********************************************************************************
C                                                                              *
C     Note:  This function checks the return value of any call to the          *
C     netcdf interface.  This subroutine should be called after every          *
C     call to any nf_* functions.                                              *
C                                                                              *
C     Author: A. Dingwell                                                      *
C                                                                              *
********************************************************************************
        implicit none
        include 'netcdf.inc'
        integer errcode

        if( errcode.ne.nf_noerr ) then
          print *, 'Error: ', nf_strerror(errcode)
          stop
        endif
      end
