      real*8 function Gauss2D(lon,lat,dt,row,col,hd)
!利用高斯基函数方法内插(lon,lat)处的函数值
!2021年4月20日,章传银
!-------------------------------------------------------------
      implicit none
	integer row,col,ii,jj,ki,kj
	real*8::dt(row,col),hd(6),fun
	real*8::qi,dm,rr,lat,lon,lat0,lon0,RAD,ae
!-------------------------------------------------------------
      RAD=datan(1.d0)/45.d0;ae=6378000
	ii=nint((lat-hd(3))/hd(6)+0.5d0)
	jj=nint((lon-hd(1))/hd(5)+0.5d0)
	qi=0.d0;dm=0.d0;Gauss2D=9999.d0
	do ki=ii-4,ii+4
	  lat0=hd(3)+(dble(ki)-0.5d0)*hd(6)
	  do kj=jj-4,jj+4
	    lon0=hd(1)+(dble(kj)-0.5d0)*hd(5)
	    rr=((lon-lon0)*dcos((lat+lat0)/2.d0*RAD))**2+(lat-lat0)**2
	    if(ki>0 .and. kj>0 .and. ki<=row .and. kj<=col)then
            if(dt(ki,kj)<9000.d0)then
             fun=dexp(-rr/hd(5)**2)
             dm=dm+dt(ki,kj)*fun;qi=qi+fun
            endif
	    endif
	  enddo
	enddo
	if(qi>1.d-18)Gauss2D=dm/qi
9100  continue
      end
!
!******************************************************************
!
      real*8 function CGrdPntD(lon,lat,dt,row,col,hd)
!利用距离反比方法内插(lon,lat)处的函数值
!2005年4月20日,章传银
!-------------------------------------------------------------
      implicit none
	integer row,col,ii,jj,ki,kj
	real*8::dt(row,col),hd(6),mdr
	real*8::qi,dm,rr,lat,lon,lat0,lon0,RAD,ae
!-------------------------------------------------------------
      RAD=datan(1.d0)/45.d0;ae=6378000
	ii=nint((lat-hd(3))/hd(6)+0.5d0)
	jj=nint((lon-hd(1))/hd(5)+0.5d0)
	qi=0.d0;dm=0.d0;CGrdPntD=9999.d0
      mdr=hd(5)*RAD*ae*1.d-3
	do ki=ii-4,ii+4
	  lat0=hd(3)+(dble(ki)-0.5d0)*hd(6)
	  do kj=jj-4,jj+4
	    lon0=hd(1)+(dble(kj)-0.5d0)*hd(5)
	    rr=((lon-lon0)*dcos((lat+lat0)/2.d0*RAD))**2+(lat-lat0)**2
	    rr=dsqrt(rr)*ae+mdr
	    if(ki>0 .and. kj>0 .and. ki<=row .and. kj<=col)then
            if(dt(ki,kj)<9000.d0)then
	        qi=qi+1.d0/rr;dm=dm+dt(ki,kj)/rr
            endif
	    endif
	  enddo
	enddo
	if(qi>1.d-16)CGrdPntD=dm/qi
9100  continue
      end
!
!******************************************************************
!
      real*8 function CGrdPntD2(lon,lat,dt,row,col,hd)
!利用距离平方反比方法内插(lon,lat)处的函数值
!2005年4月20日,章传银
!-------------------------------------------------------------
      implicit none
	integer row,col,ii,jj,ki,kj,kn
	real*8::dt(row,col),hd(6),mdr,tmp
	real*8::qi,dm,rr,lat,lon,lat0,lon0,RAD,ae
!-------------------------------------------------------------
      RAD=datan(1.d0)/45.d0;ae=6378000
	ii=nint((lat-hd(3))/hd(6)+0.5d0)
	jj=nint((lon-hd(1))/hd(5)+0.5d0)
	qi=0.d0;dm=0.d0;tmp=0.d0
      mdr=hd(5)*RAD*ae*1.d-3
	do ki=ii-4,ii+4
	  lat0=hd(3)+(dble(ki)-0.5d0)*hd(6)
	  do kj=jj-4,jj+4
	    lon0=hd(1)+(dble(kj)-0.5d0)*hd(5)
	    rr=((lon-lon0)*dcos(lat*RAD))**2+(lat-lat0)**2
	    rr=rr*ae**2+mdr**2
	    if(ki>0 .and. kj>0 .and. ki<=row .and. kj<=col)then
            if(dt(ki,kj)<9000.d0)then
	        qi=qi+1.d0/rr;dm=dm+dt(ki,kj)/rr
            endif
          endif
	  enddo
	enddo
	if(qi>1.d-16)tmp=dm/qi
      CGrdPntD2=tmp
      end
!
!******************************************************************
!
      real*8 function CShepard(lon,lat,dt,row,col,hd)
!利用Shepard方法内插(lon,lat)处的函数值
!2005年4月20日,章传银
!-------------------------------------------------------------
      implicit none
	integer row,col,ii,jj,ki,kj
	real*8::dt(row,col),hd(6),mdr,R0,pr
	real*8::qi,dm,rr,lat,lon,lat0,lon0,RAD,ae
!-------------------------------------------------------------
      RAD=datan(1.d0)/45.d0;ae=6378000
	ii=nint((lat-hd(3))/hd(6)+0.5d0)
	jj=nint((lon-hd(1))/hd(5)+0.5d0)
	qi=0.d0;dm=0.d0;CShepard=9999.d0;
      mdr=hd(5)*RAD*ae*1.d-3
	R0=4.d0*hd(5)*RAD*ae
	do ki=ii-4,ii+4
	  lat0=hd(3)+(dble(ki)-0.5d0)*hd(6)
	  do kj=jj-4,jj+4
	    lon0=hd(1)+(dble(kj)-0.5d0)*hd(5)
	    rr=(((lon-lon0)*dcos(lat*RAD))**2+(lat-lat0)**2)*RAD**2
	    rr=dsqrt(rr)*ae+mdr
	    if(ki>0 .and. kj>0 .and. ki<=row .and. kj<=col)then
	      if(rr<=R0/3.d0.and.dt(ki,kj)<9000.d0)then
	        qi=qi+(1.d0/rr)**2;dm=dm+dt(ki,kj)/rr**2
	      endif
	      if(rr>R0/3.d0.and.rr<=R0.and.dt(ki,kj)<9000.d0)then
	        pr=27.d0/4.d0/R0*(rr/R0-1.d0)**2
	        qi=qi+pr**2;dm=dm+dt(ki,kj)*pr**2
	      endif
	    endif
	  enddo
	enddo
	if(qi>1.d-12)CShepard=dm/qi
9100  continue
      end
