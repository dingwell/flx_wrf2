      subroutine writeheader()
********************************************************************************
*                                                                              *
*  Note:  This is the FLEXPART_WRF version of subroutine writeheader.          *
*                                                                              *
*  This routine produces a file header containing basic information on the     *
*  settings of the FLEXPART run.                                               *
*  The header file is essential and must be read in by any postprocessing      *
*  program before reading in the output data.                                  *
*                                                                              *
*     Author: A. Stohl                                                         *
*     7 August 2002                                                            *
*                                                                              *
*     Dec 2005, J. Fast & R. Easter -                                          *
*            Write formatted output when iouttype=1.                           *
*            Write iomode_xycoord.  Write xy coords. as lat-lon or             *
*            grid-meters depending on iomode_xycoord.                          *
*            changed names of "*lon0*" & "*lat0*" variables                    *
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
* ymet                   model y coordinate in grid-meters                     *
*                                                                              *
* xp1                    temporary x-coordinate                                *
* yp1                    temporary y-coordinate                                *
* xp2                                                                          *
* yp2                                                                          *
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
* ncsrcid                ID for netcdf sources dimension (npoints)             *
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
      real xv1(nspec),yv1(nspec),xv2(nspec),yv2(nspec)
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

      if (iouttype .eq. 0) then       ! binary
        open(unitheader,file=path(2)(1:len(2))//'header',
     +  form='unformatted',err=998)
      else if (iouttype .eq. 1) then  ! ascii
        open(unitheader,file=path(2)(1:len(2))//'header',
     +  form='formatted',err=998)
      else                            ! netcdf
        write(ncname,'(A11,I8.8,A1,I6.6,A3)')
     +  'flxout_d01_',ibdate,'_',ibtime,'.nc' ! filename
        ncret = nf_create(path(2)(1:len(2))//ncname,
     +  nf_clobber,ncid) ! open & overwrite if it exists
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
        ncret = nf_put_att_text(ncid,nf_global,
     +    'TITLE',20,version)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +    'START_DATE',nf_int,1,ibdate)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +    'START_TIME',nf_int,1,ibtime)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +    'END_DATE',nf_int,1,iedate)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +    'END_TIME',nf_int,1,ietime)
        call check_ncerror(ncret)
      endif !iouttype

C Write info on output interval, averaging time, sampling time
C FLEXPART_WRF - also the iomode_xycoord
**************************************************************

      if (iouttype .eq. 0) then       ! binary
        write(unitheader) loutstep,loutaver,loutsample,iomode_xycoord
      else if (iouttype .eq. 1) then   ! ascii
        write(unitheader,*) loutstep,loutaver,loutsample,iomode_xycoord
      else                            ! netcdf
        ncret = nf_put_att_int(ncid,nf_global,
     +  'OUTPUT_INTERVAL',nf_int,1,loutstep)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'AVERAGING_TIME',nf_int,1,loutaver)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'AVERAGE_SAMPLING',nf_int,1,loutsample)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'IOMODE_XYCOORD',nf_int,1,iomode_xycoord)
        call check_ncerror(ncret)
      endif !iouttype

C Write information on output grid setup
****************************************

      call caldate(bdate,jjjjmmdd,ihmmss)

      if (iomode_xycoord .eq. iomode_xycoord_latlon) then
         xp1 = outgrid_swlon
         yp1 = outgrid_swlat
         dxtmp = (outgrid_nelon-outgrid_swlon)/numxgrid
         dytmp = (outgrid_nelat-outgrid_swlat)/numygrid
      else
         xp1 = out_xm0
         yp1 = out_ym0
         dxtmp = dxout
         dytmp = dyout
      endif

      if (iouttype .eq. 0) then       ! binary
        write(unitheader) xp1,yp1,numxgrid,numygrid,dxtmp,dytmp
        write(unitheader) numzgrid,(outheight(i),i=1,numzgrid)
        write(unitheader) jjjjmmdd,ihmmss
      else if (iouttype .eq. 1) then  ! ascii
        write(unitheader,*) xp1,yp1,numxgrid,numygrid,dxtmp,dytmp
        write(unitheader,*) numzgrid,(outheight(i),i=1,numzgrid)
        write(unitheader,*) jjjjmmdd,ihmmss
      else                            ! netcdf
        ncret = nf_put_att_int(ncid,nf_global,
     +  'OUTLONLEFT',nf_int,1,xp1)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'OUTLATLOWER',nf_int,1,yp1)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'WEST-EAST_GRID_DIMENSION',nf_int,1,numxgrid)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'SOUTH-NORTH_GRID_DIMENSION',nf_int,1,numygrid)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
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
        ncret = nf_def_dim(ncid,  ! TIME (record)
     +  'Times',nf_unlimited,ncrecid)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncid,  ! WEST-EAST
     +  'west_east',numxgrid,nclonid)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncid,  ! NORTH-SOUTH
     +  'south_north',numygrid,nclatid)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncid,  ! BOTTOM-UP
     +  'bottom_top',numzgrid,nclvlid)
        call check_ncerror(ncret)
        write(*,*) ncrecid,nclonid,nclatid,nclvlid

        ncret = nf_def_dim(ncid,  ! SPECIES
     +  'species',nspec,ncspcid)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncid,  ! SPECIES NAME
     +  'SpeciesStrLen',10,ncstr1id)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncid,  ! AGECLASSES
     +  'ageclass',nageclass,ncageid)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncid,  ! SOURCES (NPOINTS)
     +  'sources',numpoint,ncsrcid)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncid,  ! SOURCE COMMENT
     +  'SourceStrLen',45,ncstr2id)
        call check_ncerror(ncret)

        ncret = nf_def_dim(ncid,  ! SOURCE START-END
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
        ncret = nf_def_var(ncid,
     +  'XLONG',nf_real,2,ncdimsid(1:2),nclonvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,nclonvid,descr,
     +  27,'LONGITUDE, WEST IS NEGATIVE')
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,nclonvid,units,
     +  11,'degree_east')
        call check_ncerror(ncret)

        ! Setup dimension variable XLAT
        ncret = nf_def_var(ncid,
     +  'XLAT',nf_real,2,ncdimsid(1:2),nclatvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,nclatvid,descr,
     +  27,'LATITUDE, SOUTH IS NEGATIVE')
        ncret = nf_put_att_text(ncid,nclatvid,units,
     +  12,'degree_north')
        call check_ncerror(ncret)

        ! Setup dimension variable ZTOP
        ncret = nf_def_var(ncid,
     +  'ZTOP',nf_real,1,ncdimsid(3),nclvlvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,nclvlvid,descr,
     +  32,'HEIGHT OF LEVEL (UPPER BOUNDARY)')
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,nclvlvid,units,1,'m')
        call check_ncerror(ncret)

        ! Setup dimension variable SPECIES
        ncret = nf_def_var(ncid,
     +  'SPECIES',nf_real,2,(/ncstr1id,ncspcid/),ncspcvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncspcvid,descr,
     +  15,'NAME OF SPECIES')
        call check_ncerror(ncret)

        ! Setup dimension variable AGECLASSES
        ncret = nf_def_var(ncid,
     +  'AGECLASS',nf_int,1,ncageid,ncagevid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncagevid,descr,
     +  27,'MAX AGE OF SPECIES IN CLASS')
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncagevid,units,1,'s')
        call check_ncerror(ncret)

c       ncret = nf_def_dim(ncid,nf_global,  ! AGECLASSES
c    +  'species',nageclass,ncageid)
c       call check_ncerror(ncret)
      endif !iouttype

C Write number of species, and name for each species (+extra name for depositions)
C Indicate the dimension of the fields (i.e., 1 for deposition fields, numzgrid for
C concentration fields
***********************************************************************************

      if (iouttype .eq. 0) then       ! binary
        write(unitheader) 3*nspec, numreceptor, nageclass
      else if (iouttype .eq. 1) then  ! ascii
        write(unitheader,*) 3*nspec, numreceptor, nageclass
      else                            ! netcdf
        ncret = nf_put_att_int(ncid,nf_global,
     +  'NSPEC',nf_int,1,nspec)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'NUMRECEPTOR',nf_int,1,numreceptor)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
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

C Write information on release points: total number, then for each point:
C start, end, coordinates, # of particles, name, mass
*************************************************************************

      if (iouttype .eq. 0) then       ! binary
        write(unitheader) numpoint
      else if (iouttype .eq. 1) then  ! ascii
        write(unitheader,*) numpoint
      else                            ! netcdf
        ! Define source related attributes
        ncret = nf_put_att_int(ncid,nf_global,
     +  'NUM_RELEASES',nf_int,1,numpoint)
        call check_ncerror(ncret)

        ! Define Source related variables
        ncret = nf_def_var(ncid,      ! Source Name
     +  'SourceName',nf_char,2,(/ncstr2id,ncsrcid/),ncsnvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncsnvid,descr,
     +  23,'SOURCE IDENTIER/COMMENT')
        ncret = nf_put_att_text(ncid,ncsnvid,units,1,'-')
        call check_ncerror(ncret)

        ncret = nf_def_var(ncid,      ! Source start-end time
     +  'SourceTstart_end',nf_int,2,(/ncsseid,ncsrcid/),ncstvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncstvid,descr,
     +  32,'BEGINNING/ENDING TIME OF RELEASE (SECONDS SINCE RUN START)')
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncstvid,units,
     +  1,'s')
c    +  14,'YYYYMMDDhhmmss')
        call check_ncerror(ncret)

        ncret = nf_def_var(ncid,      ! Source start-end X-coord
     +  'SourceXstart_end',nf_float,2,(/ncsseid,ncsrcid/),ncsxvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncsxvid,descr,
     +  30,'WEST/EAST BOUNDARIES OF SOURCE')
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncsxvid,units,
     +  11,'degree_east')
        call check_ncerror(ncret)

        ncret = nf_def_var(ncid,      ! Source start-end Y-coord
     +  'SourceYstart_end',nf_float,2,(/ncsseid,ncsrcid/),ncsyvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncsyvid,descr,
     +  32,'SOUTH/NORTH BOUNDARIES OF SOURCE')
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncsyvid,units,
     +  12,'degree_north')
        call check_ncerror(ncret)

        ncret = nf_def_var(ncid,      ! Source start-end Z-coord
     +  'SourceZstart_end',nf_float,2,(/ncsseid,ncsrcid/),ncszvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncszvid,descr,
     +  31,'BOTTOM/TOP BOUNDARIES OF SOURCE')
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncszvid,units,1,'m')
        call check_ncerror(ncret)

        ncret = nf_def_var(ncid,      ! Number of particles
     +  'SourceNP',nf_int,1,ncsrcid,ncspvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncspvid,descr,
     +  34,'TOTAL NUMBER OF PARTICLES RELEASED')
        ncret = nf_put_att_text(ncid,ncspvid,units,1,'-')
        call check_ncerror(ncret)

        ncret = nf_def_var(ncid,      ! Mass Emission
     +  'SourceXMass',nf_real,2,(/ncspcid,ncsrcid/),ncsmvid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncsmvid,descr,
     +  18,'Total mass emitted')
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,ncsmvid,units,2,'kg')
        call check_ncerror(ncret)

        ! Since we need to exit define mode before we can insert
        ! variable data, we will include the last file attributes and
        ! define the last variables here.

        ! OUTPUT VARIABLES (all netcdf variables must be predefined)
        ncret = nf_def_var(ncid,      ! Topography
     +  'TOPO',NF_real,2,ncdimsid(1:2),nctovid)
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,nctovid,descr,
     +  30,'TERRAIN HEIGHT ABOVE SEA LEVEL')
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,nctovid,units,1,'m')
        call check_ncerror(ncret)
        ncret = nf_put_att_text(ncid,nctovid,coord,coordxylen,coordxy)
        call check_ncerror(ncret)

        ! Check which output grids should be used
        if (ldirect .eq. 1) then  ! Forward run
          ! Deposition grids:
          ncret = nf_def_var(ncid,'DRYDEP',NF_REAL, ! Dry deposition
     +    5,(/nclonid,nclatid,ncspcid,ncageid,ncrecid/),ncddvid)
          call check_ncerror(ncret)
          ncret = nf_put_att_text(ncid,ncddvid,descr,
     +    32,'ACCUMULATED TOTAL DRY DEPOSITION')
          call check_ncerror(ncret)
          ncret = nf_put_att_text(ncid,ncddvid,units,6,'ng m-2')
          call check_ncerror(ncret)
          ncret = nf_put_att_text(ncid,ncddvid,coord,coordxylen,coordxy)
          call check_ncerror(ncret)

          ncret = nf_def_var(ncid,'WETDEP',NF_REAL, ! Wet deposition
     +    5,(/nclonid,nclatid,ncspcid,ncageid,ncrecid/),ncwdvid)
          call check_ncerror(ncret)
          ncret = nf_put_att_text(ncid,ncwdvid,descr,
     +    32,'ACCUMULATED TOTAL WET DEPOSITION')
          call check_ncerror(ncret)
          ncret = nf_put_att_text(ncid,ncwdvid,units,6,'ng m-2')
          call check_ncerror(ncret)
          ncret = nf_put_att_text(ncid,ncwdvid,coord,coordxylen,coordxy)
          call check_ncerror(ncret)

          if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
            ncret = nf_def_var(ncid,      ! Concentration
     +      'CONC',NF_REAL,6,ncdimsid,nccovid)
            call check_ncerror(ncret)
            ncret = nf_put_att_text(ncid,nccovid,descr,
     +      33,'CONCENTRATION OF AIRBORNE SPECIES')
            call check_ncerror(ncret)
            ncret = nf_put_att_text(ncid,nccovid,units,6,'ng m-3')
            call check_ncerror(ncret)
            ncret = nf_put_att_text(ncid,nccovid,coord,
     +      coordxylen,coordxy)
            call check_ncerror(ncret)
          endif
          if ((iout.eq.2).or.(iout.eq.3)) then
            ncret = nf_def_var(ncid,      ! Mixing ratio
     +      'RATIO',NF_REAL,6,ncdimsid,ncravid)
            call check_ncerror(ncret)
            ncret = nf_put_att_text(ncid,ncravid,descr,
     +      37,'MASS MIXING RATIO OF AIRBORNE SPECIES')
            call check_ncerror(ncret)
            ncret = nf_put_att_text(ncid,ncravid,units,3,'ppt')
            call check_ncerror(ncret)
            ncret = nf_put_att_text(ncid,ncravid,coord,
     +      coordxylen,coordxy)
            call check_ncerror(ncret)
          endif
        else                      ! Backward run
          if (( iout .eq. 1 ).or.( iout .eq. 5 )) then  ! Use restime?
            write(*,*)  '### NETCDF NOT SUPPORTED FOR BACKWARD RUNS ###'
            stop
            !TODO Setup residence time variable
            ! I don't think this should be included in the netCDF
            ! there will be difficulties due to the strict format
            ! (only one expandable variable)
            ! however, this information doesn't really require any
            ! advanced formating, so maybe it's not a big problem to 
            ! keep it as ascii?
          endif
          if (( iout .eq. 4 ).or.( iout .eq. 5 )) then ! Use traj
            write(*,*)  '### NETCDF NOT SUPPORTED FOR BACKWARD RUNS ###'
            stop
            !TODO Setup trajectory var
          endif
        endif !ldirect


        ! MODEL SWITCHES (NETCDF ONLY, ascii+bin is further down)
        ncret = nf_put_att_int(ncid,nf_global,
     +  'DISPERSION_METHOD',nf_int,1,method)
        call check_ncerror(ncret)
        
        ncret = nf_put_att_int(ncid,nf_global,
     +  'SUBGRID_TOPOGRAPHY',nf_int,1,lsubgrid)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'CONVECTION_PARAMETERIZATION',nf_int,1,lconvection)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'SUBGRID_TOPOGRAPHY',nf_int,1,lsubgrid)
        call check_ncerror(ncret)

        ncret = nf_put_att_int(ncid,nf_global,
     +  'NUM_RELEASES',nf_int,1,numpoint)
        call check_ncerror(ncret)

        ! Exit netcdf define mode, enter data mode
        ncret = nf_enddef(ncid)
        call check_ncerror(ncret)

        ! Dimension variables (netcdf)
        ! Z-height
        ncret = nf_put_var_real(ncid,nclvlvid,outheight)
        call check_ncerror(ncret)
        ! X,Y-Lon,Lat
        call ll_to_xymeter_wrf(outgrid_swlon,outgrid_swlat,xsw,ysw)
        call ll_to_xymeter_wrf(outgrid_nelon,outgrid_nelat,xne,yne)
        do jy=1,numygrid
        do ix=1,numxgrid
          tmpx=xsw+(xne-xsw)*float(ix-1)/float(numxgrid-1)
          tmpy=ysw+(yne-ysw)*float(jy-1)/float(numygrid-1)
          call xymeter_to_ll_wrf(tmpx,tmpy,tmplon,tmplat)
          ncret = nf_put_vara_real(ncid,nclonvid,
     +    (/ix,jy/),(/1,1/),tmplon)
          call check_ncerror(ncret)
          ncret = nf_put_vara_real(ncid,nclatvid,
     +    (/ix,jy/),(/1,1/),tmplat)
          call check_ncerror(ncret)
        enddo
        enddo
      endif !iouttype

      do 13 i=1,numpoint  ! (loop is used for binary and ascii only)

C New method
        if (iomode_xycoord .eq. iomode_xycoord_latlon) then ! latlon
          xv1(i)=releases_swlon(i)    ! use ll as is
          yv1(i)=releases_swlat(i)    !
          xv2(i)=releases_nelon(i)    !
          yv2(i)=releases_nelat(i)    !
        else                                                ! metres
          xv1(i)=xpoint1(i)*dx+xmet0  ! easting (metres)
          yv1(i)=ypoint1(i)*dx+ymet0  ! northing(metres)
          xv2(i)=xpoint2(i)*dx+xmet0  !
          yv2(i)=ypoint2(i)*dx+ymet0  !
        endif
C OLD METHOD FOR RELEASE COORDINATES KEPT JUST IN CASE, REMOVE LATER?
C       if (iomode_xycoord .eq. iomode_xycoord_latlon) then ! latlon
C          xp1=releases_swlon(i)    ! use ll as is
C          yp1=releases_swlat(i)    !
C          xp2=releases_nelon(i)    !
C          yp2=releases_nelat(i)    !
C       else                                                ! metres
C          xp1=xpoint1(i)*dx+xmet0  ! easting (metres)
C          yp1=ypoint1(i)*dy+ymet0  ! northing(metres)
C          xp2=xpoint2(i)*dx+xmet0  !
C          yp2=ypoint2(i)*dy+ymet0  !
C       endif

        if (iouttype .eq. 0) then     ! binary
          write(unitheader) ireleasestart(i),ireleaseend(i),kindz(i)
          write(unitheader) 
     +     xv1(i),yv1(i),xv2(i),yv2(i),zpoint1(i),zpoint2(i)
          write(unitheader) npart(i),1
          write(unitheader) compoint(i)
          do j=1,nspec
            write(unitheader) xmass(i,j)
            write(unitheader) xmass(i,j)
            write(unitheader) xmass(i,j)
          enddo
        else if (iouttype .eq. 1) then! ascii
          write(unitheader,*) ireleasestart(i),ireleaseend(i),kindz(i)
          write(unitheader,*) 
     +      xv1(i),yv1(i),xv2(i),yv2(i),zpoint1(i),zpoint2(i)
          write(unitheader,*) npart(i),1
          write(unitheader,*) compoint(i)
          do j=1,nspec
            write(unitheader,*) xmass(i,j)
            write(unitheader,*) xmass(i,j)
            write(unitheader,*) xmass(i,j)
          enddo
        else                          ! netcdf
          ncret = nf_put_vara_int(ncid,ncstvid,   ! SourceTstart_end
     +    (/1,i/),(/2,1/),(/ireleasestart(i),ireleaseend(i)/))
          call check_ncerror(ncret)

          ncret = nf_put_vara_real(ncid,ncsxvid,  ! SourceXstart_end
     +    (/1,i/),(/2,1/),(/xv1(i),xv2(i)/))
          call check_ncerror(ncret)

          ncret = nf_put_vara_real(ncid,ncsyvid,  ! SourceYstart_end
     +    (/1,i/),(/2,1/),(/yv1(i),yv2(i)/))
          call check_ncerror(ncret)

          ncret = nf_put_vara_real(ncid,ncszvid,  ! SourceZstart_end
     +    (/1,i/),(/2,1/),(/zpoint1(i),zpoint2(i)/))
          call check_ncerror(ncret)

          ncret = nf_put_vara_real(ncid,ncsmvid,  ! SourceXMass
     +    (/1,i/),(/nspec,1/),xmass(i,1:nspec))
          call check_ncerror(ncret)

          ncret = nf_put_vara_int(ncid,ncspvid,  ! SourceNP
     +    i,1,npart(i))
          call check_ncerror(ncret)

          !Source Name/Comment
          j=1 ! Find the length of each source comment/name
          do while( j.lt.45.and.compoint(i)(j+1:j+1).ne." ")
            j=j+1
          end do
            ncret = nf_put_vara_text(ncid,ncsnvid,   ! write to file
     +      (/1,i/),(/j,1/),compoint(i)(1:j))
            call check_ncerror(ncret)
        endif

13    continue  ! do..numpoint

C Write information on some model switches
******************************************

      if (iouttype .eq. 0) then
        write(unitheader) method,lsubgrid,lconvection
      else if (iouttype .eq. 1) then
        write(unitheader,*) method,lsubgrid,lconvection
      endif

C Write age class information
*****************************

      if (iouttype .eq. 0) then
        write(unitheader) nageclass,(lage(i),i=1,nageclass)
      else if (iouttype .eq. 1) then
        write(unitheader,*) nageclass,(lage(i),i=1,nageclass)
      else
        ncret = nf_put_var_int(ncid,ncagevid,lage(1:nageclass))
        call check_ncerror(ncret)
      endif


C Write topography to output file
*********************************

      do 30 ix=0,numxgrid-1
        if (iouttype .eq. 0) then
          write(unitheader) (oroout(ix,jy),jy=0,numygrid-1)
        else if (iouttype .eq. 1) then
          write(unitheader,*) (oroout(ix,jy),jy=0,numygrid-1)
        else 
          ncret = nf_put_vara_real(ncid,nctovid,
     +    (/ix+1,1/),(/1,numygrid/),oroout(ix,0:numygrid-1))
          call check_ncerror(ncret)
        endif

      if(iouttype.eq.2) then          ! netcdf
        ncret = nf_sync(ncid)   ! Save changes to file
        call check_ncerror(ncret)
      endif
30    continue

      close(unitheader)
      return


998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
      write(*,*) ' #### '//path(2)(1:len(2))//'header'//' #### '
      write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
      write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
      write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
      stop

      end
