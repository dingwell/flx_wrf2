      subroutine writeheader_nest()
********************************************************************************
*                                                                              *
*  Note:  This is the FLEXPART_WRF version of subroutine writeheder_nest.      *
*                                                                              *
*  This routine produces a file header containing basic information on the     *
*  settings of the FLEXPART run.                                               *
*  The header file is essential and must be read in by any postprocessing      *
*  program before reading in the output data.                                  *
*                                                                              *
*     Author: A. Stohl                                                         *
*     7 August 2002                                                            *
*                                                                              *
*     Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables     *
*                                                                              *
*     Nov 2012, A. Dingwell -                                                  *
*            Setup NetCDF object and write to file when iouttype=2             *
*            The same file will be used by concoutput.f and is referenced by   *
*            the global variable ncid                                          *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
* xmet                   model x coordinate in grid-meters                     *
* ymet                   model x coordinate in grid-meters                     *
*                                                                              *
* xv1                    temporary array of x-coordinates                      *
* yv1                    temporary array of y-coordinates                      *
* xv2                                                                          *
* yv2                                                                          *
*                                                                              *
* ncret                  return value when using the netcdf Fortran interface  *
*                                                                              *
* nclvlid                ID for netcdf output vertical dimension               *
* nclonid                ID for netcdf output longitude dimension              *
* nclatid                ID for netcdf output latitude dimension               *
* ncspcid                ID for netcdf output species dimension                *
* ncageid                ID for netcdf output ageclass dimension               *
* ncrecid                ID for netcdf output time (record) dimension          *
*                                                                              *
* ncstr1id               ID for netcdf species names dimension (string)        *
*                                                                              *
* nclvlvid               ID for netcdf output vertical dimension variable      *
* nclonvid               ID for netcdf output longitude dimension variable     *
* nclatvid               ID for netcdf output latitude dimension variable      *
* ncspcvid               ID for netcdf output species dimension variable       *
* ncagevid               ID for netcdf output ageclass dimension variable      *
*                                                                              *
* ncsrcid                ID for netcdf sources dimension (numpoint)             *
* ncstr2id               ID for netcdf source names/comments dimension (string)*
* ncsseid                ID for netcdf source start_end dimension (2)          *
*                                                                              *
* ncsnvid                ID for netcdf output source name variable             *
* ncsmvid                ID for netcdf output source XMass variable            *
* ncspvid                ID for netcdf output source particles variable        *
* ncstvid                ID for netcdf output source time-span variable        *
* ncsxvid                ID for netcdf output source x-span variable           *
* ncsyvid                ID for netcdf output source y-span variable           *
* ncszvid                ID for netcdf output source z-span variable           *
*                                                                              *
* nctovid                ID for netcdf output topography variable              *
* ncdimsid               Array with netcdf dimension ids (lon,lat,lvl,time)    *
* ncname                 filename of output file (only for netCDF)             *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'
      include 'netcdf.inc'
      
      integer jjjjmmdd,ihmmss,i,ix,jy,j

      real dxtmp,dytmp
      real xp1,yp1,xp2,yp2  ! use vectors instead:
      real xv1(numpoint),yv1(numpoint),xv2(numpoint),yv2(numpoint)
      real xsw,xne,ysw,yne,tmpx,tmpy,tmplon,tmplat  ! get ll of outgrid

      ! NETCDF SPECIFIC VARIABLES
      integer nclvlid,nclonid,nclatid,ncrecid,ncspcid,ncageid ! outgrid dimensions
      integer ncsrcid,ncsseid                     ! source dimensions
      integer ncsnvid,ncsmvid,ncspvid             ! source strength def. vars
      integer ncstvid,ncsxvid,ncsyvid,ncszvid     ! source box defining vars
      integer ncstr1id,ncstr2id                   ! string length dims
      integer nclvlvid,nclonvid,nclatvid,ncspcvid,ncagevid    ! outgrid dim-variables
      integer ncdimsid(6)
      character descr*11,units*5,ncname*29,coord*11,coordxy*10
      integer coordxylen
      integer nctovid ! local netcdf variables
      
      descr = 'description'
      units = 'units'
      coord = 'coordinates'
      coordxy = 'XLONG XLAT'
      coordxylen = 10


*************************
C Open header output file
*************************

      if(iouttype.eq.0) then          ! binary
        open(unitheader,file=path(2)(1:len(2))//'header_nest',
     +  form='unformatted',err=998)
      else if (iouttype .eq. 1) then  ! ascii
        open(unitheader,file=path(2)(1:len(2))//'header',
     +  form='formatted',err=998)
      else                            ! netcdf
        write(ncname,'(A8,I2.2,A1,I8.8,A1,I6.6,A3)')
     +  'flxout_d',i+1,'_',ibdate,'_',ibtime,'.nc' ! filename
        ncret = nf_create(path(2)(1:len(2))//ncname,
     +  nf_clobber,ncidn) ! open & overwrite if it exists
        call check_ncerror(ncret)
      endif


C Write the header information
******************************

      if(iouttype.eq.0.or.iouttype.eq.1) then ! binary or ascii
        if (ldirect.eq.1) then  ! forward run
          if (iouttype .eq. 0) then     ! binary
            write(unitheader) ibdate,ibtime,'FLEXPWRF V5.0'
          else if (iouttype.eq.1) then  ! ascii
            write(unitheader,*) ibdate,ibtime
            write(unitheader,*) 'FLEXPWRF V5.0'
          endif !iouttype
        else                    ! backward run
          if (iouttype .eq. 0) then     ! binary
            write(unitheader) iedate,ietime,'FLEXPWRF V5.0'
          else if (iouttype .eq. 1) then! ascii
            write(unitheader,*) iedate,ietime
            write(unitheader,*) 'FLEXPWRF V5.0'
          endif !iouttype
        endif !ldirect
      else if(iouttype.eq.2) then ! netcdf
        ncret = nf_put_att_text(ncidn,nf_global,
     +    'TITLE',20,version)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncidn,nf_global,
     +    'START_DATE',nf_int,1,ibdate)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncidn,nf_global,
     +    'START_TIME',nf_int,1,ibtime)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncidn,nf_global,
     +    'END_DATE',nf_int,1,iedate)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncidn,nf_global,
     +    'END_TIME',nf_int,1,ietime)
        call check_ncerror(ncret)
      endif !iouttype

C Write info on output interval, averaging time, sampling time
**************************************************************

      if (iouttype .eq. 0) then       ! binary
        write(unitheader) loutstep,loutaver,loutsample,iomode_xycoord
      else if (iouttype .eq. 1) then   ! ascii
        write(unitheader,*) loutstep,loutaver,loutsample,iomode_xycoord
      else                            ! netcdf
        ncret = nf_put_att_int(ncidn,nf_global,
     +  'OUTPUT_INTERVAL',nf_int,1,loutstep)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncidn,nf_global,
     +  'AVERAGING_TIME',nf_int,1,loutaver)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncidn,nf_global,
     +  'AVERAGE_SAMPLING',nf_int,1,loutsample)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncidn,nf_global,
     +  'IOMODE_XYCOORD',nf_int,1,iomode_xycoord)
        call check_ncerror(ncret)
      endif !iouttype

C Write information on output grid setup
****************************************

      call caldate(bdate,jjjjmmdd,ihmmss)

      if (iouttype.eq.0) then ! binary
        write(unitheader) out_xm0n,out_ym0n,
     +  numxgridn,numygridn,
     +  dxoutn,dyoutn
        write(unitheader) numzgrid,(outheight(i),i=1,numzgrid)
        write(unitheader) jjjjmmdd,ihmmss
      else if (iouttype.eq.1) then ! ascii
        write(unitheader,*) out_xm0n,out_ym0n,
     +  numxgridn,numygridn,
     +  dxoutn,dyoutn
        write(unitheader,*) numzgrid,(outheight(j),j=1,numzgrid)
        write(unitheader,*) jjjjmmdd,ihmmss
      else                    ! netcdf
        ncret = nf_put_att_int(ncid,nf_global,
     +  'OUTLONLEFT',nf_int,1,xp1)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'OUTLATLOWER',nf_int,1,yp1)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncidn,nf_global,
     +  'WEST-EAST_GRID_DIMENSION',nf_int,1,numxgridn)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncidn,nf_global,
     +  'SOUTH-NORTH_GRID_DIMENSION',nf_int,1,numygridn)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncidn,nf_global,
     +  'BOTTOM-TOP_GRID_DIMENSION',nf_int,1,numzgrid)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'DX',nf_int,1,dxtmp)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'DY',nf_int,1,dytmp)
        call check_ncerror(ncret)

        ! Not sure why these two are used, I include them until I know:
        ncret = nf_put_att_int(ncid,nf_global,
     +  'jjjjmmdd',nf_int,1,jjjjmmdd)
        ncret = nf_put_att_int(ncid,nf_global,
     +  'ihmmss',nf_int,1,ihmmss)

        ! Setup netcdf dimensions
        ncret = nf_def_dim(ncidn,  ! TIME (record)
     +  'Times',nf_unlimited,ncrecid)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncidn,  ! WEST-EAST
     +  'west_east',numxgrid,nclonid)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncidn,  ! NORTH-SOUTH
     +  'south_north',numygrid,nclatid)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncidn,  ! BOTTOM-UP
     +  'bottom_top',numzgrid,nclvlid)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncidn,  ! SPECIES
     +  'species',nspec,ncspcid)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncidn,  ! SPECIES NAME
     +  'SpeciesStrLen',10,ncstr1id)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncidn,  ! AGECLASSES
     +  'ageclass',nageclass,ncageid)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncidn,  ! SOURCES (NPOINTS)
     +  'sources',numpoint,ncsrcid)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncidn,  ! SOURCE COMMENT
     +  'SourceStrLen',45,ncstr2id)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncidn,  ! SOURCE START-END
     +  'SourceStartEnd',2,ncsseid)
        call check_ncerror(ncret)

        ! Setup array for outgrid dimensions
        ! These are the dimensions used for the main 
        ! concentration grids
        ncdimsid(1) = nclonid ! X
        ncdimsid(2) = nclatid ! Y
        ncdimsid(3) = nclvlid ! Z
        ncdimsid(4) = ncspcid ! species
        ncdimsid(5) = ncageid ! ageclass
        ncdimsid(6) = ncrecid ! t

        ! Setup dimension variable XLONG
        ncret = nf_def_var(ncidn,
     +  'XLONG',nf_real,2,ncdimsid(1:2),nclonvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncidn,nclonvid,descr,
     +  27,'LONGITUDE, WEST IS NEGATIVE')
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncidn,nclonvid,units,
     +  11,'degree_east')
        call check_ncerror(ncret)

        ! Setup dimension variable XLAT
        ncret = nf_def_var(ncidn,
     +  'XLAT',nf_real,2,ncdimsid(1:2),nclatvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncidn,nclatvid,descr,
     +  27,'LATITUDE, SOUTH IS NEGATIVE')
        ncret = nf_put_att_text(ncidn,nclatvid,units,
     +  12,'degree_north')
        call check_ncerror(ncret)

        ! Setup dimension variable ZTOP
        ncret = nf_def_var(ncidn,
     +  'ZTOP',nf_real,1,ncdimsid(3),nclvlvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncidn,nclvlvid,descr,
     +  32,'HEIGHT OF LEVEL (UPPER BOUNDARY)')
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncidn,nclvlvid,units,1,'m')
        call check_ncerror(ncret)

        ! Setup dimension variable SPECIES
        ncret = nf_def_var(ncidn,
     +  'SPECIES',nf_char,2,(/ncstr1id,ncspcid/),ncspcvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncidn,ncspcvid,descr,
     +  15,'NAME OF SPECIES')
        call check_ncerror(ncret)

        ! Setup dimension variable AGECLASSES
        ncret = nf_def_var(ncidn,
     +  'AGECLASS',nf_int,1,ncageid,ncagevid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncidn,ncagevid,descr,
     +  27,'MAX AGE OF SPECIES IN CLASS')
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncidn,ncagevid,units,1,'s')
        call check_ncerror(ncret)

c       ncret = nf_def_dim(ncid,nf_global,  ! AGECLASSES
c      +'species',nageclass,ncageid)
c       call check_ncerror(ncret)
      endif !iouttype


C Write number of species, and name for each species (+extra name for depositions)
C Indicate the dimension of the fields (i.e., 1 for deposition fields, numzgrid for
C concentration fields
***********************************************************************************

      if (iouttype .eq. 0) then       ! binary
        write(unitheader) 3*nspec
      else if (iouttype .eq. 1) then  ! ascii
        write(unitheader,*) 3*nspec
      else                            ! netcdf
        ncret = nf_put_att_int(ncidn,nf_global,
     +  'NSPEC',nf_int,1,nspec)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'NUMRECEPTOR',nf_int,1,numreceptor)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncidn,nf_global,
     +  'NAGECLASS',nf_int,1,nageclass)
        call check_ncerror(ncret)
      endif

      do 12 i=1,nspec
         if (iouttype .eq. 0) then      ! binary
           write(unitheader) 1,'WD_'//species(i)(1:7)
           write(unitheader) 1,'DD_'//species(i)(1:7)
           write(unitheader) numzgrid,species(i)
         else if (iouttype .eq. 1) then ! ascii
           write(unitheader,*) 1
           write(unitheader,*) 'WD_'//species(i)(1:7)
           write(unitheader,*) 1
           write(unitheader,*) 'DD_'//species(i)(1:7)
           write(unitheader,*) numzgrid
           write(unitheader,*) species(i)
         endif
12    continue

C We will not write information on release points, this information is 
C already available in the main output file
*************************************************************************

C Write information on some model switches
******************************************

      if (iouttype .eq. 0) then      ! binary
        write(unitheader) method,lsubgrid,lconvection
      else if (iouttype .eq. 1) then ! ascii or netcdf
        write(unitheader,*) method,lsubgrid,lconvection
      endif
C AD: sorry for not finishing this, but it was starting to use up too
C much of my time.  This is as far I got. In order to make nests work for
C netcdf we must create new variable pointers for the nested fields,
C otherwise they will overwrite the ones used by the mother domain.
C When writeheader_nest.f is finally working, additional changes are
C required in concoutput_nest.f similar to my previous modifications in
C concoutput.f. I might take care of this later, but I have no project
C planned right now which would benefit from it, so it might take some
C time.

C Write age class information
*****************************

      write(unitheader) nageclass,(lage(i),i=1,nageclass)


C Write topography to output file
*********************************

      do 30 ix=0,numxgridn-1
30      write(unitheader) (orooutn(ix,jy),jy=0,numygridn-1)
      close(unitheader)

      return


998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
      write(*,*) ' #### '//path(2)(1:len(2))//'header'//' #### '
      write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
      write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
      write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
      stop

      end
