      subroutine subdag(nm,n,low,igh,h,wr,wi,ierr,norm,k,its,en,na,enm2,
     2 l,s,t,retVal,x,y,w,itn) 
      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
      integer retVal
      double precision h(nm,n),wr(n),wi(n)
      double precision p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
      logical notlas

c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
      do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = dabs(h(l-1,l-1)) + dabs(h(l,l))
         if (s .eq. 0.0d0) s = norm
         tst1 = s
         tst2 = tst1 + dabs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
  100 return 
      end
