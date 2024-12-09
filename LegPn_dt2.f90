       subroutine LegPn_dt2(pn,dp1,dp2,n,t)
      !Legendre function and its first and second derivatives to θ. 
      !计算Pn(t)及其对θ一、二阶导数t=cosθ
	implicit none
	integer::n,k
	real*8::pn(n),dp1(n),dp2(n),t,u
!---------------------------------------------------------------------------
      u=dsqrt(1-t**2)
      pn(1)=t;pn(2)=1.5d0*t**2-0.5d0
      dp1(1)=-u;dp1(2)=-3.d0*t*u
      dp2(1)=-t;dp2(2)=3.d0*(1.d0-2.d0*t**2)
      do k=3,n
        pn(k)=dble(2*k-1)/dble(k)*t*pn(k-1)-dble(k-1)*pn(k-2)/dble(k)
        dp1(k)=dble(2*k-1)*(t*dp1(k-1)-u*pn(k-1))-dble(k-1)*dp1(k-2)
        dp1(k)=dp1(k)/dble(k)
        dp2(k)=(t*dp2(k-1)-2.d0*u*dp1(k-1)-t*pn(k-1))*dble(2*k-1)
        dp2(k)=(dp2(k)-dble(k-1)*dp2(k-2))/dble(k)
      enddo
      end
