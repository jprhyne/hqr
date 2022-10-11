      subroutine tmphqr(nm,n,low,igh,h,wr,wi,ierr,norm,k,its,en,na,enm2,
     2 ll,l,s) 
      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
      double precision h(nm,n),wr(n),wi(n)
      double precision p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
      logical notlas
