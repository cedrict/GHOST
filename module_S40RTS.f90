module s40rts_tomography

integer,parameter :: MXLH=40
integer,parameter :: MXLENY=(MXLH+1)**2
integer,parameter :: MXDEP=21
integer,parameter :: MXSPL=MXDEP+3
integer,parameter :: MXMSZ=MXSPL*MXLENY
integer,parameter :: MAXP=1024
real(8),parameter ::  xi_s40rts=.25

real wk1(MXLENY),wk2(MXLENY),wk3(MXLENY)
real d0(MXLENY)
real spl(MXDEP)
real x(MXMSZ)
real dum,xd

real, parameter :: rcmb=3480.d3
real, parameter :: rmoho=6346.d3
real, parameter :: rearth=6371.d3

integer nbeg,nend,ipp,ind,ind1
integer lmax,natd,ndmx,nsmx,ndmn

end module
