      subroutine frmsft(nm,n,low,igh,h,wr,wi,ierr,norm,k,its,en,na,enm2,
     2 l,s,t,retVal,x,y,w,itn,m) 
      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
      integer retVal
      double precision h(nm,n),wr(n),wi(n)
      double precision p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
      logical notlas

      x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
c     .......... form exceptional shift ..........
      t = t + x
c
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
c
      s = dabs(h(en,na)) + dabs(h(na,enm2))
      x = 0.75d0 * s
      y = x
      w = -0.4375d0 * s * s
      go to 1001
  270 retVal = 1
      go to 1001
  280 retVal = 2
      go to 1001
  130 go to 1001
 1000 ierr = en
 1001 return
      end
