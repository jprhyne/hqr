      subroutine hqr_subdiagsearch(nm,n,low,h,norm,l,s) 
      implicit none
      integer ll,low,en,l
      double precision h(nm,n)
      double precision s,norm
c---------------------------------------------------------------
c       Variables Modified: ll,l,s
c       Variables only read: low,en,h,norm
c       Variables used locally only: tst1,tst2
c---------------------------------------------------------------

c     local variables
      double precision tst1,tst2

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
