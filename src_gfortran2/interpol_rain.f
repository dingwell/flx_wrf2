      subroutine interpol_rain(yy1,yy2,yy3,nxmax,nymax,nzmax,nx,
     +ny,memind,xt,yt,level,itime1,itime2,itime,yint1,yint2,yint3)
C                               i   i   i    i    i     i   i 
C     i    i    i  i    i     i      i      i     o     o     o
*****************************************************************************
*                                                                           *
*  Interpolation of meteorological fields on 2-d model layers.              *
*  In horizontal direction bilinear interpolation interpolation is used.    *
*  Temporally a linear interpolation is used.                               *
*  Three fields are interpolated at the same time.                          *
*                                                                           *
*  This is a special version of levlininterpol to save CPU time.            *
*                                                                           *
*  1 first time                                                             *
*  2 second time                                                            *
*                                                                           *
*                                                                           *
*     Author: A. Stohl                                                      *
*                                                                           *
*     30 August 1996                                                        *
*                                                                           *
*     19 December 2012: Adam Dingwell                                       *
*       yy1 and yy2 are only interpolated in space, while yy3 is still      *
*       interpolated in both space and time.  This is due to the way        *
*       precipitation is read from the WRF model, there should be no        *
*       interpolation between time steps for lsprec or convprec!            *
*       The old parts of this file are kept as comments in case a better    *
*       way of getting precipitation data is included in a future version.  *
*                                                                           *
*****************************************************************************
*                                                                           *
* Variables:                                                                *
*                                                                           *
* dt1,dt2              time differences between fields and current position *
* dz1,dz2              z distance between levels and current position       *
* height(nzmax)        heights of the model levels                          *
* indexh               help variable                                        *
* indz                 the level closest to the current trajectory position *
* indzh                help variable                                        *
* itime                current time                                         *
* itime1               time of the first wind field                         *
* itime2               time of the second wind field                        *
* ix,jy                x,y coordinates of lower left subgrid point          *
* level                level at which interpolation shall be done           *
* memind(3)            points to the places of the wind fields              *
* nx,ny                actual field dimensions in x,y and z direction       *
* nxmax,nymax,nzmax    maximum field dimensions in x,y and z direction      *
* xt                   current x coordinate                                 *
* yint                 the final interpolated value                         *
* yt                   current y coordinate                                 *
* yy(0:nxmax,0:nymax,nzmax,3) meteorological field used for interpolation   *
* zt                   current z coordinate                                 *
*                                                                           *
*****************************************************************************

      implicit none

      integer nx,ny,nxmax,nymax,nzmax,memind(2),m,ix,jy,ixp,jyp
      integer itime,itime1,itime2,level,indexh
      real yy1(0:nxmax-1,0:nymax-1,nzmax,2)
      real yy2(0:nxmax-1,0:nymax-1,nzmax,2)
      real yy3(0:nxmax-1,0:nymax-1,nzmax,2)
      real ddx,ddy,rddx,rddy,dt1,dt2,dt,y1(2),y2(2),y3(2)
      real xt,yt,yint1,yint2,yint3,p1,p2,p3,p4



C If point at border of grid -> small displacement into grid
************************************************************

      if (xt.ge.float(nx-1)) xt=float(nx-1)-0.00001
      if (yt.ge.float(ny-1)) yt=float(ny-1)-0.00001

C AD: Do the same for domain lower boundary:
      if (xt.le.0.) xt=0.00001
      if (yt.le.0.) yt=0.00001



***********************************************************************
C 1.) Bilinear horizontal interpolation
C This has to be done separately for 2 fields (Temporal)
********************************************************

C Determine the lower left corner and its distance to the current position
************************************************************************** 

      ix=int(xt)
      jy=int(yt)
      ixp=ix+1
      jyp=jy+1
      ddx=xt-float(ix)
      ddy=yt-float(jy)
      rddx=1.-ddx
      rddy=1.-ddy
      p1=rddx*rddy
      p2=ddx*rddy
      p3=rddx*ddy
      p4=ddx*ddy
     

C Loop over 2 time steps
************************
 
C AD: Comment out the following 8 lines to enable
C     temporal interpolation of precipitation
      yint1 = p1*yy1(ix ,jy ,level,memind(2))
     +      + p2*yy1(ixp,jy ,level,memind(2))
     +      + p3*yy1(ix ,jyp,level,memind(2))
     +      + p4*yy1(ixp,jyp,level,memind(2))
      yint2 = p1*yy2(ix ,jy ,level,memind(2))
     +      + p2*yy2(ixp,jy ,level,memind(2))
     +      + p3*yy2(ix ,jyp,level,memind(2))
     +      + p4*yy2(ixp,jyp,level,memind(2))
C /AD

      do 10 m=1,2
        indexh=memind(m)
        
C AD: Uncomment the following 8 lines to enable
C     temporal interpolation of precipitation
c       y1(m)=p1*yy1(ix ,jy ,level,indexh)
c    +      + p2*yy1(ixp,jy ,level,indexh)
c    +      + p3*yy1(ix ,jyp,level,indexh)
c    +      + p4*yy1(ixp,jyp,level,indexh)
c       y2(m)=p1*yy2(ix ,jy ,level,indexh)
c    +      + p2*yy2(ixp,jy ,level,indexh)
c    +      + p3*yy2(ix ,jyp,level,indexh)
c    +      + p4*yy2(ixp,jyp,level,indexh)
C /AD
10      y3(m)=p1*yy3(ix ,jy ,level,indexh)
     +      + p2*yy3(ixp,jy ,level,indexh)
     +      + p3*yy3(ix ,jyp,level,indexh)
     +      + p4*yy3(ixp,jyp,level,indexh)


*************************************
C 2.) Temporal interpolation (linear)
*************************************

      dt1=float(itime-itime1)
      dt2=float(itime2-itime)
      dt=dt1+dt2

C AD: Uncomment the following 2 lines to enable
C     temporal interpolation of precipitation
c     yint1=(y1(1)*dt2+y1(2)*dt1)/dt
c     yint2=(y2(1)*dt2+y2(2)*dt1)/dt
C /AD
      yint3=(y3(1)*dt2+y3(2)*dt1)/dt

      if(yint1.lt.0 .or. yint2.lt.0) then
        write(*,*) 'interpol_rain: NEGATIVE PRECIP! yint1,yint2:',
     +  yint1,yint2
        write(*,*) 'yy1:',yy1(ix ,jy ,level,memind(2))
     +                   ,yy1(ixp,jy ,level,memind(2))
     +                   ,yy1(ix ,jyp,level,memind(2))
     +                   ,yy1(ixp,jyp,level,memind(2))
        write(*,*) 'p1,p2,p3,p4:',p1,p2,p3,p4
        write(*,*) 'xt,yt:',xt,yt
        write(*,*) 'ddx,ddy,rddx,rddy:',ddx,ddy,rddx,rddy
      end if

      end
