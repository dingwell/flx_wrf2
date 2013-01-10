      subroutine wetdepo(itime,ltsample,loutnext)
C                          i      i        i
********************************************************************************
*                                                                              *
* Calculation of wet deposition using the concept of scavenging coefficients.  *
* For lack of detailed information, washout and rainout are jointly treated.   *
* It is assumed that precipitation does not occur uniformly within the whole   *
* grid cell, but that only a fraction of the grid cell experiences rainfall.   *
* This fraction is parameterized from total cloud cover and rates of large     *
* scale and convective precipitation.                                          *
*                                                                              *
*    Author: A. Stohl                                                          *
*                                                                              *
*    1 December 1996                                                           *
*                                                                              *
* Correction by Petra Seibert, Sept 2002:                                      *
* use centred precipitation data for integration                               *
* Code may not be correct for decay of deposition!                             *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* cc [0-1]           total cloud cover                                         *
* convp [mm/h]       convective precipitation rate                             *
* fraction [0-1]     fraction of grid, for which precipitation occurs          *
* ix,jy              indices of output grid cell for each particle             *
* itime [s]          actual simulation time [s]                                *
* jpart              particle index                                            *
* ldeltat [s]        interval since radioactive decay was computed             *
* lfr, cfr           area fraction covered by precipitation for large scale    *
*                    and convective precipitation (dependent on prec. rate)    *
* loutnext [s]       time for which gridded deposition is next output          *
* loutstep [s]       interval at which gridded deposition is output            *
* lsp [mm/h]         large scale precipitation rate                            *
* ltsample [s]       interval over which mass is deposited                     *
* prec [mm/h]        precipitation rate in subgrid, where precipitation occurs *
* wetdeposit         mass that is wet deposited                                *
* wetgrid            accumulated deposited mass on output grid                 *
* wetscav            scavenging coefficient                                    *
*                                                                              *
* Constants:                                                                   *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer jpart,itime,ltsample,loutnext,ldeltat,i,j,k,ix,jy
      integer ngrid,itage,nage,hz,il,interp_time,n,clouds_v
      real S_i, act_temp, cl, cle ! in cloud scavenging
      real clouds_h ! cloud height for the specific grid point
      real xtn,ytn,lsp,convp,cc,fraction,prec,wetscav
      real wetdeposit(maxspec),lfr(5),cfr(5),restmass,smallnum
      save lfr,cfr,smallnum

      data lfr/0.5,0.65,0.8,0.9,0.95/
      data cfr/0.4,0.55,0.7,0.8,0.9/
      data smallnum /1.e-38/ ! should be smallest number that can be handled

C Compute interval since radioactive decay of deposited mass was computed
*************************************************************************

      if (itime.le.loutnext) then
        ldeltat=itime-(loutnext-loutstep)
      else                                  ! first half of next interval
        ldeltat=itime-loutnext
      endif


C Loop over all particles
*************************

      do 20 jpart=1,numpart
cAF        if (itra1(jpart).eq.-999999999) goto 20
         if(ldirect.eq.1)then
            if (itra1(jpart).gt.itime) goto 20
         else
            if (itra1(jpart).lt.itime) goto 20
         endif   
C Determine age class of the particle
        itage=abs(itra1(jpart)-itramem(jpart))
        do 32 nage=1,nageclass
          if (itage.lt.lage(nage)) goto 33
32        continue
33      continue


C Determine which nesting level to be used
******************************************

        ngrid=0
        do 22 j=numbnests,1,-1
          if ((xtra1(jpart).gt.xln(j)).and.(xtra1(jpart).lt.xrn(j)).and.
     +    (ytra1(jpart).gt.yln(j)).and.(ytra1(jpart).lt.yrn(j))) then
            ngrid=j
            goto 23
          endif
22        continue
23      continue


C Determine nested grid coordinates
***********************************

        if (ngrid.gt.0) then
          xtn=(xtra1(jpart)-xln(ngrid))*xresoln(ngrid)
          ytn=(ytra1(jpart)-yln(ngrid))*yresoln(ngrid)
          ix=int(xtn)
          jy=int(ytn)
        else
          ix=int(xtra1(jpart))
          jy=int(ytra1(jpart))
        endif


C Interpolate large scale precipitation, convective precipitation and
C total cloud cover
C Note that interpolated time refers to itime-0.5*ltsample [PS]
*********************************************************************
        interp_time=nint(itime-0.5*ltsample)

        if (ngrid.eq.0) then
          call interpol_rain(lsprec,convprec,tcc,nxmax,nymax,
     +    1,nx,ny,memind,sngl(xtra1(jpart)),sngl(ytra1(jpart)),1,
     +    memtime(1),memtime(2),interp_time,lsp,convp,cc)
        else
          call interpol_rain_nests(lsprecn,convprecn,tccn,
     +    nxmaxn,nymaxn,1,maxnests,ngrid,nxn,nyn,memind,xtn,ytn,1,
     +    memtime(1),memtime(2),interp_time,lsp,convp,cc)
        endif
C       write(*,*) 'ngrid:',ngrid
        if(lsp.lt.0.or.convp.lt.0) then
          write(*,*) 'wetdepo (1): NEGATIVE PRECIP! lsp,convp:',
     +    lsp,convp
        end if

        if ((lsp.lt.0.01).and.(convp.lt.0.01)) goto 20

C AD: adopted from FLEXPART-9.02:
C     Get the level where the particle is in
        do il=2,nz
          if (height(il).gt.ztra1(jpart)) then
            hz=il-1   ! level found
            goto 26   ! exit loop
          endif
        end do
26      continue

        n=memind(2)
        if(abs(memtime(1)-interp_time).lt.abs(memtime(2)-interp_time))
     +    n=memind(1)

C AD: adopted from FLEXPART-9.02:
C     if there is no precipitation or the particle is above the clouds
C     no scavenging is done
        
        if(ngrid.eq.0) then
          clouds_v=clouds(ix,jy,hz,n)
          clouds_h=cloudsh(ix,jy,1,n)
        else
          clouds_v=cloudsn(ix,jy,hz,n,ngrid)
          clouds_h=cloudsnh(ix,jy,1,n,ngrid)
        endif
        if (clouds_v.le.1) goto 20  ! ABOVE CLOUD
C       write (*,*) 'there is scavenging'
C /AD
        

C 1) Parameterization of the the area fraction of the grid cell where the
C    precipitation occurs: the absolute limit is the total cloud cover, but
C    for low precipitation rates, an even smaller fraction of the grid cell
C    is used. Large scale precipitation occurs over larger areas than
C    convective precipitation.
***************************************************************************

        if (lsp.gt.20.) then
          i=5
        else if (lsp.gt.8.) then
          i=4
        else if (lsp.gt.3.) then
          i=3
        else if (lsp.gt.1.) then
          i=2
        else
          i=1
        endif

        if (convp.gt.20.) then
          j=5
        else if (convp.gt.8.) then
          j=4
        else if (convp.gt.3.) then
          j=3
        else if (convp.gt.1.) then
          j=2
        else
          j=1
        endif

        if(lsp.lt.0.or.convp.lt.0) then
          write(*,*) 'wetdepo (2): NEGATIVE PRECIP! lsp,convp:',
     +    lsp,convp
        end if

        fraction=max(0.05,cc*(lsp*lfr(i)+convp*cfr(j))/(lsp+convp))

C 2) Computation of precipitation rate in sub-grid cell
*******************************************************

        prec=(lsp+convp)/fraction

        if(lsp.lt.0.or.convp.lt.0) then
          write(*,*) 'wetdepo (3): NEGATIVE PRECIP! lsp,convp:',
     +    lsp,convp
        end if

C 3) Computation of scavenging coefficients for all species
C    Computation of wet deposition
C AD: FLEXPART-9.02 code included to deal with in/below
C     cloud scavenging and disable above cloud scavenging.
***********************************************************

        do 10 k=1,nspec                               ! loop over species

        if(lsp.lt.0.or.convp.lt.0) then
          write(*,*) 'wetdepo (4): NEGATIVE PRECIP! lsp,convp:',
     +    lsp,convp
        end if

          wetdeposit(k)=0.                            ! AD: default = no scav
          if (weta(k).gt.0.) then ! species affected by wetdep?
            if (clouds_v.ge.4) then ! BELOW CLOUD?
              wetscav=weta(k)*prec**wetb(k)           ! scavenging coeff.
              ! for aerosols and not highliy soluble substances
              ! weta=5E-6
            else  ! IN CLOUD? (above cloud has already been ruled out)
              if (ngrid.gt.0) then
                act_temp=ttn(ix,jy,hz,n,ngrid)
              else
                act_temp=tt(ix,jy,hz,n)
              endif
              cl=2E-7*prec**0.36
              if(lsp.lt.0.or.convp.lt.0) then
                write(*,*) 'lsp,convp,fraction:',lsp,convp,fraction
                write(*,*) 'NEGATIVE PRECIP'
                stop
              end if
              if(dquer(k).gt.0) then ! PARTICLE
                !write(*,*) 'is particle'
                S_i=0.9/cl
              else                    ! GAS
                !write(*,*) 'is gas'
                cle=(1-cl)/(henry(k)*(r_air/3500.)*act_temp)+cl
                S_i=1/cle
              endif
              wetscav=S_i*prec/3.6E6/clouds_h
c             write(*,*) 'S_i:',S_i
c             write(*,*) 'in. wetscav:',
c    +        wetscav,cle,cl,act_temp,prec,clouds_h
            endif
            wetdeposit(k)=xmass1(jpart,k)*
     +        (1.-exp(-wetscav*abs(ltsample)))*fraction  ! wet deposition
C           new particle mass:
            restmass = xmass1(jpart,k)-wetdeposit(k)
            if (restmass .gt. smallnum) then
              xmass1(jpart,k)=restmass
            else
              xmass1(jpart,k)=0.
            endif


C Correct deposited mass to the last time step when radioactive decay of 
C gridded deposited mass was calculated
************************************************************************
            if (decay(k).gt.0.) then
              wetdeposit(k)=wetdeposit(k)*exp(abs(ldeltat)*decay(k))
            endif
          else
            wetdeposit(k)=0.
          endif
10        continue

C AD: Followed Sabine Eckhard's example: only call wetdepokernel for
C forward runs.
C Add the wet deposition to accumulated amount on output grid and nested output grid
************************************************************************************

        if (ldirect.eq.1) then
          call wetdepokernel(nclass(jpart),wetdeposit,sngl(xtra1(jpart))
     +    ,sngl(ytra1(jpart)),nage)
          if (nested_output.eq.1) call wetdepokernel_nest(nclass(jpart),
     +      wetdeposit,sngl(xtra1(jpart)),sngl(ytra1(jpart)),nage)
        endif

20      continue

      end
