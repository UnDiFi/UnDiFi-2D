      subroutine pspltm(nrow,ncol,mode,ja,ia,title,ptitle,size,munt,
     *     nlines,lines,iunt)
c-----------------------------------------------------------------------
      integer nrow,ncol,ptitle,mode,iunt, ja(*), ia(*), lines(nlines) 
      real size
      character title*(*), munt*2 
c----------------------------------------------------------------------- 
c PSPLTM - PostScript PLoTer of a (sparse) Matrix
c This version by loris renggli (renggli@masg1.epfl.ch), Dec 1991
c and Youcef Saad 
c------
c Loris RENGGLI, Swiss Federal Institute of Technology, Math. Dept
c CH-1015 Lausanne (Switzerland)  -- e-mail:  renggli@masg1.epfl.ch
c Modified by Youcef Saad -- June 24, 1992 to add a few features:
c separation lines + acceptance of MSR format.
c-----------------------------------------------------------------------
c input arguments description :
c
c nrow   = number of rows in matrix
c
c ncol   = number of columns in matrix 
c
c mode   = integer indicating whether the matrix is stored in 
c           CSR mode (mode=0) or CSC mode (mode=1) or MSR mode (mode=2) 
c
c ja     = column indices of nonzero elements when matrix is
c          stored rowise. Row indices if stores column-wise.
c ia     = integer array of containing the pointers to the 
c          beginning of the columns in arrays a, ja.
c
c title  = character*(*). a title of arbitrary length to be printed 
c          as a caption to the figure. Can be a blank character if no
c          caption is desired.
c
c ptitle = position of title; 0 under the drawing, else above
c
c size   = size of the drawing  
c
c munt   = units used for size : 'cm' or 'in'
c
c nlines = number of separation lines to draw for showing a partionning
c          of the matrix. enter zero if no partition lines are wanted.
c
c lines  = integer array of length nlines containing the coordinates of 
c          the desired partition lines . The partitioning is symmetric: 
c          a horizontal line across the matrix will be drawn in 
c          between rows lines(i) and lines(i)+1 for i=1, 2, ..., nlines
c          an a vertical line will be similarly drawn between columns
c          lines(i) and lines(i)+1 for i=1,2,...,nlines 
c
c iunt   = logical unit number where to write the matrix into.
c----------------------------------------------------------------------- 
c additional note: use of 'cm' assumes european format for paper size
c (21cm wide) and use of 'in' assumes american format (8.5in wide).
c The correct centering of the figure depends on the proper choice. Y.S.
c-----------------------------------------------------------------------
c external 
      integer LENSTR
      external LENSTR
c local variables ---------------------------------------------------
      integer n,nr,nc,maxdim,istart,ilast,ii,k,ltit
      real lrmrgn,botmrgn,xtit,ytit,ytitof,fnstit,siz
      real xl,xr, yb,yt, scfct,u2dot,frlw,delt,paperx,conv,xx,yy
      logical square 
c change square to .true. if you prefer a square frame around
c a rectangular matrix
      data haf /0.5/, zero/0.0/, conv/2.54/,square/.false./
c-----------------------------------------------------------------------
      siz = size
      nr = nrow
      nc = ncol
      n = nc
      if (mode .eq. 0) n = nr
c      nnz = ia(n+1) - ia(1) 
      maxdim = max(nrow, ncol)
      m = 1 + maxdim
      nc = nc+1
      nr = nr+1
c
c units (cm or in) to dot conversion factor and paper size
c 
      if (munt.eq.'cm' .or. munt.eq.'CM') then
         u2dot = 72.0/conv
        paperx = 21.0
      else
        u2dot = 72.0
        paperx = 8.5*conv
        siz = siz*conv
      end if
c
c left and right margins (drawing is centered)
c 
      lrmrgn = (paperx-siz)/2.0
c
c bottom margin : 2 cm
c
      botmrgn = 2.0
c scaling factor
      scfct = siz*u2dot/m
c matrix frame line witdh
      frlw = 0.25
c font size for title (cm)
      fnstit = 0.5
      ltit = LENSTR(title)
c position of title : centered horizontally
c                     at 1.0 cm vertically over the drawing
      ytitof = 1.0
      xtit = paperx/2.0
      ytit = botmrgn+siz*nr/m + ytitof
c almost exact bounding box
      xl = lrmrgn*u2dot - scfct*frlw/2
      xr = (lrmrgn+siz)*u2dot + scfct*frlw/2
      yb = botmrgn*u2dot - scfct*frlw/2
      yt = (botmrgn+siz*nr/m)*u2dot + scfct*frlw/2
      if (ltit.gt.0) then
        yt = yt + (ytitof+fnstit*0.70)*u2dot
      end if
c add some room to bounding box
      delt = 10.0
      xl = xl-delt
      xr = xr+delt
      yb = yb-delt
      yt = yt+delt
c
c correction for title under the drawing
      if (ptitle.eq.0 .and. ltit.gt.0) then
      ytit = botmrgn + fnstit*0.3
      botmrgn = botmrgn + ytitof + fnstit*0.7
      end if
c begin of output
c
      write(iunt,10) '%!'
      write(iunt,10) '%%Creator: PSPLTM routine'
      write(iunt,12) '%%BoundingBox:',xl,yb,xr,yt
      write(iunt,10) '%%EndComments'
      write(iunt,10) '/cm {72 mul 2.54 div} def'
      write(iunt,10) '/mc {72 div 2.54 mul} def'
      write(iunt,10) '/pnum { 72 div 2.54 mul 20 string'
      write(iunt,10) 'cvs print ( ) print} def'
      write(iunt,10)
     1  '/Cshow {dup stringwidth pop -2 div 0 rmoveto show} def'
c
c we leave margins etc. in cm so it is easy to modify them if
c needed by editing the output file
      write(iunt,10) 'gsave'
      if (ltit.gt.0) then
      write(iunt,*) '/Helvetica findfont ',fnstit,
     &             ' cm scalefont setfont '
      write(iunt,*) xtit,' cm ',ytit,' cm moveto '
      write(iunt,'(3A)') '(',title(1:ltit),') Cshow'
      end if
      write(iunt,*) lrmrgn,' cm ',botmrgn,' cm translate'
      write(iunt,*) siz,' cm ',m,' div dup scale '
c------- 
c draw a frame around the matrix
      write(iunt,*) frlw,' setlinewidth'
      write(iunt,10) 'newpath'
      write(iunt,11) 0, 0, ' moveto'
      if (square) then
      write(iunt,11) m,0,' lineto'
      write(iunt,11) m, m, ' lineto'
      write(iunt,11) 0,m,' lineto'
      else
      write(iunt,11) nc,0,' lineto'
      write(iunt,11) nc,nr,' lineto'
      write(iunt,11) 0,nr,' lineto'
      end if
      write(iunt,10) 'closepath stroke'
c
c     drawing the separation lines 
c 
      write(iunt,*)  ' 0.2 setlinewidth'
      do 22 kol=1, nlines 
         isep = lines(kol) 
c
c     horizontal lines 
c
         yy =  real(nrow-isep) + haf 
         xx = real(ncol+1) 
         write(iunt,13) zero, yy, ' moveto '
         write(iunt,13)  xx, yy, ' lineto stroke '
c
c vertical lines 
c
         xx = real(isep) + haf 
         yy = real(nrow+1)  
         write(iunt,13) xx, zero,' moveto '
         write(iunt,13) xx, yy, ' lineto stroke '             
 22     continue
c 
c----------- plotting loop ---------------------------------------------
c
      write(iunt,10) '1 1 translate'
      write(iunt,10) '0.8 setlinewidth'
      write(iunt,10) '/p {moveto 0 -.40 rmoveto '
      write(iunt,10) '           0  .80 rlineto stroke} def'
c     
      do 1 ii=1, n
        istart = ia(ii)
        ilast  = ia(ii+1)-1 
        if (mode .eq. 1) then
          do 2 k=istart, ilast
            write(iunt,11) ii-1, nrow-ja(k), ' p'
 2        continue 
        else
          do 3 k=istart, ilast
            write(iunt,11) ja(k)-1, nrow-ii, ' p'
 3        continue          
c add diagonal element if MSR mode.
          if (mode .eq. 2) 
     *         write(iunt,11) ii-1, nrow-ii, ' p' 
c
        endif
 1    continue
c-----------------------------------------------------------------------
      write(iunt,10) 'showpage'
      return
c
 10   format (A)
 11   format (2(I6,1x),A)
 12   format (A,4(1x,F9.2))
 13   format (2(F9.2,1x),A)
c-----------------------------------------------------------------------
      end
