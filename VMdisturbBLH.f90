      subroutine VMdisturbBLH(BLH,gra,hgt,nlat,nlon,hd,dr,vm,GRS)
      !按严密球面积分公式计算广义V_M扰动重力积分
      !dr-积分半径m
!-------------------------------------------------------------
      implicit none
	integer::i,j,nlat,nlon,i0,j0,ni,nj
	real*8::dr,gra(nlat,nlon),hgt(nlat,nlon),vm(2),rst(2)
	real*8::GRS(6),hd(6),mr,pi,RAD,ds,mdr,tt,rr,r0,r1,r2
	real*8::BLH(3),XYZ(3),rln(3),BLH0(3),XYZ0(3),BLH1(3),XYZ1(3),rln1(3)
	real*8 CGrdPntD2,HotineSgn,L0,L1,SK,gr,NFD(5)
	real*8 rlat,rlon,rlat1,rlon1,sinf,sin2f,cosa,sina
!-----------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0;mr=3.6d-2/RAD;BLH0=BLH
      BLH0(3)=CGrdPntD2(BLH(2),BLH(1),hgt,nlat,nlon,hd)!计算点正下方地面点位置
      call BLH_RLAT(GRS,BLH,rln);rr=rln(1);call BLH_XYZ(GRS,BLH,XYZ)
      call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      call BLH_XYZ(GRS,BLH0,XYZ0)
      r0=dsqrt(XYZ0(1)**2+XYZ0(2)**2+XYZ0(3)**2)
      mdr=r0*hd(5)*RAD*dcos(rln(2)*RAD)/2.d0 !奇异点判断
      ni=nint(dr/r0/RAD/hd(6)+1.d0) !积分半径dr对应的地面格网数
      nj=nint(dr/r0/RAD/hd(5)/dcos(rln(2)*RAD)+1.d0)
	i0=nint((BLH(1)-hd(3))/hd(6)+0.5d0)
	j0=nint((BLH(2)-hd(1))/hd(5)+0.5d0)!计算点所在的地面格网i0,j0
      rlat=rln(2)*RAD;rlon=rln(3)*RAD;vm=0.d0
	do i=i0-ni,i0+ni
	  if(i<1.or.i>nlat)goto 9100
        BLH1(1)=hd(3)+(real(i)-0.5d0)*hd(6)
	  do j=j0-nj,j0+nj
	    if(j<1.or.j>nlon)goto 9101
	    BLH1(2)=hd(1)+(real(j)-0.5d0)*hd(5)
          BLH1(3)=hgt(i,j);call BLH_XYZ(GRS,BLH1,XYZ1)
          L0=dsqrt((XYZ1(1)-XYZ0(1))**2+(XYZ1(2)-XYZ0(2))**2+(XYZ1(3)-XYZ0(3))**2)
          if(L0>dr)goto 9101
          L1=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          if(L0<mdr)then!计算奇异积分
             call VMdturbSgn(BLH,gra,hgt,nlat,nlon,hd,rst,i,j,4,GRS)
             vm(1)=vm(1)+rst(1);vm(2)=vm(2)+rst(2)
             goto 9101 
          endif
          tt=1.d0-2.d0*(L0/r1/2.d0)**2
          rlat1=rln1(2)*RAD;rlon1=rln1(3)*RAD
          ds=hd(5)*hd(6)*RAD**2*dcos(rlat1)*r1**2
          sin2f=dsqrt((1.d0-tt)/2.d0);sinf=2.d0*sin2f*dsqrt((1.d0+tt)/2.d0)
	    cosa=(dcos(rlat)*dsin(rlat1)-dsin(rlat)*dcos(rlat1)*dcos(rlon1-rlon))/sinf
	    sina=dcos(rlat1)*dsin(rlon1-rlon)/sinf
          SK=-2.d0*rr*r1/L1**3+3.d0*r1/rr**2-(L1-rr)/L1/(rr+L1-r1*tt)
          SK=SK+1.d0/rr/(1.d0-tt);SK=SK*sinf
	    vm(1)=vm(1)+gra(i,j)*SK*cosa*ds/4.d0/pi/gr*mr/rr
	    vm(2)=vm(2)+gra(i,j)*SK*sina*ds/4.d0/pi/gr*mr/rr
9101      continue
	  enddo
9100    continue
	enddo
9002	return
      end
!--------------------------------------------------------------------------------
      subroutine VMdturbSgn(BLH,gra,hgt,nlat,nlon,hd,vm,i0,j0,m,GRS)
      !细化核函数，计算BLH点的广义V_M扰动重力奇异积分
      !m-核函数细化为m*m
      !i0,j0-奇异点格网位置
!-------------------------------------------------------------
      implicit none
	integer::m,i,j,nlat,nlon,i0,j0,ni,nj
	real*8::dr,gra(nlat,nlon),hgt(nlat,nlon),vm(2)
	real*8::GRS(6),hd(6),pi,RAD,ds,mdr,tt,rr,r0,r1,r2,rv,mr
	real*8::BLH(3),XYZ(3),rln(3),BLH0(3),XYZ0(3),rln0(3),BLH1(3),XYZ1(3),rln1(3)
	real*8 CGrdPntD2,L0,L1,SK,lon,lat,dg,gr,NFD(5)
	real*8 rlat,rlon,rlat1,rlon1,sinf,sin2f,cosa,sina
!-----------------------------------------------------------------
      BLH0=BLH;rv=hd(5)/dble(m);pi=datan(1.d0)*4.d0;RAD=pi/180.d0;mr=3.6d-2/RAD
      lat=hd(3)+real(i0-1)*hd(6);lon=hd(1)+real(j0-1)*hd(5)!格网左下角经纬度
      BLH0(3)=CGrdPntD2(BLH(2),BLH(1),hgt,nlat,nlon,hd)!计算点正下方地面点位置
      call BLH_XYZ(GRS,BLH,XYZ);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      rr=dsqrt(XYZ(1)**2+XYZ(2)**2+XYZ(3)**2)
      call BLH_XYZ(GRS,BLH0,XYZ0);call BLH_RLAT(GRS,BLH0,rln0);r0=rln0(1)
      mdr=r0*rv*RAD*dcos(rln0(2)*RAD)/4.d0 !奇异点判断
      rlat=rln0(2)*RAD;rlon=rln0(3)*RAD;vm=0.d0
	do i=1,m
        BLH1(1)=lat+(real(i)-0.5d0)*rv
	  do j=1,m
	    BLH1(2)=lon+(real(j)-0.5d0)*rv
          BLH1(3)=CGrdPntD2(BLH1(2),BLH1(1),hgt,nlat,nlon,hd)
          dg=CGrdPntD2(BLH1(2),BLH1(1),gra,nlat,nlon,hd)
          call BLH_XYZ(GRS,BLH1,XYZ1)
          L0=dsqrt((XYZ1(1)-XYZ0(1))**2+(XYZ1(2)-XYZ0(2))**2+(XYZ1(3)-XYZ0(3))**2)
          L1=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          if(L1<mdr)L1=mdr
          tt=dsqrt(1.d0-(L0/r1)**2)
          rlat1=rln1(2)*RAD;rlon1=rln1(3)*RAD
          ds=rv**2*RAD**2*dcos(rlat1)*r1**2
	    cosa=(dcos(rlat)*dsin(rlat1)-dsin(rlat)*dcos(rlat1)*dcos(rlon1-rlon))/sinf
	    sina=dcos(rlat1)*dsin(rlon1-rlon)/sinf
          SK=-2.d0*rr*r1/L1**3+3.d0*r1/rr**2-(L1-rr)/L1/(rr+L1-r1*tt)
          SK=SK+1.d0/rr/(1.d0-tt);SK=SK*sinf
	    vm(1)=vm(1)+dg*SK*cosa*ds/4.d0/pi/gr*mr/rr
	    vm(2)=vm(2)+dg*SK*sina*ds/4.d0/pi/gr*mr/rr
9101      continue
	  enddo
9100    continue
	enddo
9002	return
      end
