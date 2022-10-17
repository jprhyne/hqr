      subroutine hqr_dblSubDiag(h,enm2,
     2 l,s,x,y,w,p,q,r,zz) 
      implicit none
      integer i,l,m,mm,mp2,enm2
      double precision h(nm,n),wr(n),wi(n)
      double precision p,q,r,s,w,x,y,zz,tst1,tst2
c-------------------------------------------------------------------
c       Variables modified: m,mm,zz,r,s,p,q
c       Variables only read: l,enm2,x,y,w,h
c       Variables only locally used: tst1,tst2,mp2,i
c-------------------------------------------------------------------
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = dabs(p) + dabs(q) + dabs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = dabs(p)*(dabs(h(m-1,m-1)) + dabs(zz) + dabs(h(m+1,m+1)))
         tst2 = tst1 + dabs(h(m,m-1))*(dabs(q) + dabs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue
c
  150 mp2 = m + 2
c
      do 160 i = mp2, en
         h(i,i-2) = 0.0d0
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0d0
  160 continue
         return
         end
