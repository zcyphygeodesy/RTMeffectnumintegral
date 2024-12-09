      subroutine Stat1d0(dt,nn,rst)
      implicit none
	integer::nn,i
      real*8::dt(nn),rst(4),pv,std,maxt,mint
!---------------------------------------------------------------------
	pv=0.d0;std=0.d0;maxt=-9.d28;mint=9.d28
      do i=1,nn
        pv=pv+dt(i)
	  if(maxt<dt(i))maxt=dt(i)
        if(mint>dt(i))mint=dt(i)
      enddo
      pv=pv/real(nn)
      do i=1,nn
        std=std+dt(i)**2
      enddo
      std=sqrt(std/real(nn))
      rst(1)=pv;rst(2)=std;rst(3)=mint;rst(4)=maxt
      end
!******************************************************************
      subroutine Stat1d(dt,nn,rst)
      implicit none
	integer::nn,i,kk
      real*8::dt(nn),rst(4),pv,std,maxt,mint
!---------------------------------------------------------------------
	pv=0.d0;std=0.d0;maxt=-9.d28;mint=9.d28
      kk=0
      do i=1,nn
        if(dt(i)>9900.d0)goto 1001
        kk=kk+1;pv=pv+dt(i)
	  if(maxt<dt(i))maxt=dt(i)
        if(mint>dt(i))mint=dt(i)
1001    continue
      enddo
      pv=pv/dble(kk)
      do i=1,nn
        if(dt(i)>9900.d0)goto 1002
        std=std+(dt(i)-pv)**2
1002    continue
      enddo
      std=dsqrt(std/dble(kk))
      rst(1)=pv;rst(2)=std;rst(3)=mint;rst(4)=maxt
      end
!******************************************************************
      subroutine Stat1dnonzero(dt,nn,rst)
      implicit none
	integer::nn,i,kk
      real*8::dt(nn),rst(4),pv,std,maxt,mint
!---------------------------------------------------------------------
	pv=0.d0;std=0.d0;maxt=-9.d28;mint=9.d28
      kk=0
      do i=1,nn
        if(dt(i)>9900.d0.or.dabs(dt(i))<1.d-20)goto 1001
        kk=kk+1;pv=pv+dt(i)
	  if(maxt<dt(i))maxt=dt(i)
        if(mint>dt(i))mint=dt(i)
1001    continue
      enddo
      pv=pv/dble(kk)
      do i=1,nn
        if(dt(i)>9900.d0)goto 1002
        std=std+(dt(i)-pv)**2
1002    continue
      enddo
      std=dsqrt(std/dble(kk))
      rst(1)=pv;rst(2)=std;rst(3)=mint;rst(4)=maxt
      end
!******************************************************************
      subroutine Stat1d1(dt,nn,rst)
      implicit none
	integer::nn,i,kk
      real*8::dt(nn),rst(4),pv,std,maxt,mint
!---------------------------------------------------------------------
	pv=0.d0;std=0.d0;maxt=-9.d28;mint=9.d28
      kk=0
      do i=1,nn
        kk=kk+1;pv=pv+dt(i)
	  if(maxt<dt(i))maxt=dt(i)
        if(mint>dt(i))mint=dt(i)
      enddo
      pv=pv/dble(kk)
      do i=1,nn
        std=std+(dt(i)-pv)**2
      enddo
      std=dsqrt(std/dble(kk))
      rst(1)=pv;rst(2)=std;rst(3)=mint;rst(4)=maxt
      end
!******************************************************************
      subroutine StatGrid(grid,row,col,rst)
!-------------------------------------------------------------
      implicit none
	integer::row,col,i,j,kk
	real*8::grid(row,col),rst(4),maxt,mint,pv,err,fr
!-----------------------------------------------------------------------
	pv=0.d0;err=0.d0;maxt=-9.d28;mint=9.d28
      fr=dble(row*col);kk=0
      do i=1,row
        do j=1,col
          if(grid(i,j)>9900.d0)goto 1001
          kk=kk+1;pv=pv+grid(i,j)/fr
	    if(maxt<grid(i,j))maxt=grid(i,j)
          if(mint>grid(i,j))mint=grid(i,j)
1001      continue
        enddo
      enddo
      pv=pv/dble(kk)*fr
      do i=1,row
        do j=1,col
          if(grid(i,j)>9900.d0)goto 1002
          err=err+(grid(i,j)-pv)**2/dble(kk)
1002      continue
        enddo
      enddo
      err=dsqrt(err)
	rst(1)=pv; rst(2)=err; rst(3)=mint;rst(4)=maxt
	return
      end
!******************************************************************
      subroutine Statvect(v1,v2,row,col,rst)
!-------------------------------------------------------------
      implicit none
	integer::row,col,i,j,kk
	real*8::v1(row,col),v2(row,col),rst(4),maxt,mint,pv,err,fr
!-----------------------------------------------------------------------
	pv=0.d0;err=0.d0;maxt=-9.d28;mint=9.d28
      fr=dble(row*col);kk=0
      do i=1,row
        do j=1,col
          if(v1(i,j)>9900.d0)goto 1001
          kk=kk+1;pv=pv+v1(i,j)/fr
	    if(maxt<v1(i,j))maxt=v1(i,j)
          if(mint>v1(i,j))mint=v1(i,j)
1001      continue
        enddo
      enddo
      do i=1,row
        do j=1,col
          if(v2(i,j)>9900.d0)goto 2001
          kk=kk+1;pv=pv+v2(i,j)/fr
	    if(maxt<v2(i,j))maxt=v2(i,j)
          if(mint>v2(i,j))mint=v2(i,j)
2001      continue
        enddo
      enddo
      pv=pv/dble(kk)*fr/2.d0
      do i=1,row
        do j=1,col
          if(v1(i,j)>9900.d0)goto 1002
          err=err+(v1(i,j)-pv)**2/dble(kk)
1002      continue
        enddo
      enddo
      do i=1,row
        do j=1,col
          if(v2(i,j)>9900.d0)goto 2002
          err=err+(v2(i,j)-pv)**2/dble(kk)
2002      continue
        enddo
      enddo
      err=dsqrt(err)
	rst(1)=pv; rst(2)=err; rst(3)=mint;rst(4)=maxt
	return
      end
