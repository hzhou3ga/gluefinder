PROGRAM  GlueFinder
! pv<0.05 examines the number of identical residues and the number of cases
PARAMETER(NMAXP=2000)
PARAMETER(NMAX=2000)

PARAMETER(NCASES=98000)
PARAMETER(NCASES_TAR=1000)
!
CHARACTER*1 :: ADIR(5)
CHARACTER*3 :: AA(0:19)!LIG_LIB(500,200000)
CHARACTER*3 :: MET(2000),LIG_TYPE,MET_D(5000)
CHARACTER*5 :: TARG,NAME5
CHARACTER*12 :: NAME12,NAMEPDB(NCASES_TAR),NAME12T
CHARACTER*12 :: NAMEPARENT(NCASES)
CHARACTER*255 :: ALIGNMENT,ALIGNDIR,image
DATA ADIR/'A','B','C','D','E'/

LOGICAL*1 :: TASS(2000)

DIMENSION :: X(6)
DIMENSION :: CA_PAR(3,NMAX,NCASES),CA(3,NMAX),XHET(3,300,NCASES)
DIMENSION :: XROT(3,300)
DIMENSION :: prec(0:6,6,2),rec(0:6,6,2),anpv(0:6,6,2)

DIMENSION :: prec_met(5000),rec_met(5000),anpv_met(5000)
INTEGER	  :: NI_D(5000),NEX_D(5000),ITEM(5000)
REAL*8	:: PV_D(5000)

INTEGER :: MM1A(NMAXP),MM2A(NMAXP)
INTEGER :: JCODE(NMAXP)!,JLIG(1000)
INTEGER :: ICODE(NMAX,NCASES),NEX(0:20,2000),NHIST(0:10,0:10)
INTEGER :: NLIG(2000)
INTEGER :: NM(100)
INTEGER :: NH_PAR(NCASES)!,NTOT(1000)!,NGLIG(1000)
INTEGER :: NRES(NCASES),MET_TYPE(NCASES)!LIG_TOT(1500)
INTEGER :: NI_MET(2000)

REAL*8 PV,PV_MET(0:6,300,2000),PV_SEL(2)
DOUBLE PRECISION :: R_1(3,NMAX),R_2(3,NMAX),DRMS,W(NMAX)
DOUBLE PRECISION :: U(3,3),T(3)

DATA AA/ 'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE','TYR','TRP'/
!	*********************************************************************************************


alignment='alignmentbottom_E3_'
image='/local/images/E3ligases-20250716a/'

read(5,*)targ  
open(unit=9,action='write',file='stat_cases_v4_precision.'//targ)
	
open(unit=1,action='read',file='readjusted_precision_vs_pv_v5')
read(1,*)NDATA_p
do idum=1,NDATA_p
 read(1,*)p,r,apv,NMD,NI,NEX_I,idd,pv
        222     format(3(1f7.3,1x),4(i4,1x),1pd9.3)
                IF(PV<=0.005D0)then
			ip=1
                   	pv_SEL(1)=0.005
                else
                    ip=2
                        pv_SEL(2)=0.05
		end if
                prec(ni,nex_i,idd)=p
 	        rec(ni,nex_I,idd)=r
                anpv(ni,nex_I,idd)=apv
end do
close(1)


w=1.

open(unit=1,action='read',file='LIST.templates_ligands_all.unique.sorted')
    ! there are the ligands in the ligand template library
rewind(1)         
read(1,*)NMET 
do im=1,NMET
	read(1,*)MET(im)
end do
close(1)
    
!	-----------------------
!      Read in template structures and their associated ligands:
  
!	open(unit=20,action='read',file='LIST.'//trim(adjustl(targ))//'_1700') ! list of metabolites being screened
    open(unit=20,action='read',file='LIST.templates_ligands_with_pockets')
	rewind(20)
	read(20,*)npar
	do ik=1,npar
		read(20,*)nameparent(ik)!,id1,id2,id3
		name12=nameparent(ik)
		LIG_TYPE=name12(7:9)
		Do IM=1,NMET
			if(LIG_TYPE==MET(IM))then
				MET_TYPE(IK)=IM
				exit
			end if
		end do
 	    open(unit=54, action='read',file=trim(adjustl(image))//name12(2:3)//'/'//trim(adjustl(name12))//'.hetcord')
		rewind(54)
        read(54,*)nh_par(ik)
        nh_par(ik)=min(nh_par(ik),300)
	    do id=1,nh_par(ik)
                 read(54,*)idr,(xhet(j,id,ik),j=1,3)!,het_par(id,ik)
		end do
		close(54)

		open(unit=45,action='read',file=trim(adjustl(image))//name12(2:3)//'/'//name12(1:5)//'.cacbfiles')
			rewind(45)
			read(45,*)nres(ik)
			do i=1,Nres(ik)
				read(45,*)id,(ca_par(j,i,ik),j=1,3),(X(J),j=1,3),icode(i,ik)
			end do
			close(45)
	end do
	close(20)
	
	open(unit=43,action='read',file='LIST.'//trim(adjustl(targ))//'_2000')
	rewind(43)
	read(43,*)NPDB
	do IPDB=1,NPDB
		read(43,*)NAMEPDB(IPDB)
		read(43,*)
	end do
	close(43)
	write(6,*)'NPDB=',NPDB

! what are the ligands that bind to the given native structure?
NTOT_PREC=0
do 2999 IK=1,NPDB   
		name12T=namepdb(ik)
  open(unit=29,action='write',file='top_'//trim(adjustl(targ))//'/'//trim(adjustl(name12t))//'.top_5')

		write(6,*)IK,name12T
		PV_MET=20.
		NEX=0
		NLIG=0
		NI_MET=0
        TASS=.false.
		open(unit=45,action='read',file='cacbfiles_'//trim(adjustl(targ))//'/'//trim(adjustl(name12T))//'.cacbfiles')
		rewind(45)
		read(45,*)name5,MRES
		do i=1,MRES
			read(45,*)jd,(ca(j,i),j=1,3),(x(j),j=1,3),jcode(i)! CA coordinates, CB coordiantes, residue identity
		end do
		close(45)

		
		DO IM=1,NMET
			NLIG(IM)=0 ! now generalized for multiple ligands which ones bind to the target structure?
		end do
		DO 201 ILOOP=1,5
			aligndir=trim(adjustl(alignment))//trim(adjustl(targ))//'_'//adir(iloop)//'/'	
			open(unit=591,action='read',file=trim(adjustl(aligndir))//trim(adjustl(namePDB(IK)))//'.aln1') ! output from apoc algorithm
			rewind(591)
			read(591,*,END=201,ERR=201)ntem
!	write(6,*)'ntem=',ntem
		do 200 lk=1,ntem
     		read(591,592,END=201,ERR=201)name12,naln,nr2,nr1,ps,ncv,tm,pidtem,psim,pv
592     	format(a12,1x,3(i4,1x),1f7.3,1x,i4,1x,3(1f7.3,1x),1pd9.3)
      		read(591,593)(mm1a(i),mm2a(i),i=1,naln)
593			format(20(i5,1x,i5,1x))
			if(ps<0.2) go to 200
            if(pv>0.05D0)go to 200
			JKP=0
			do ipar=1,NPAR
				if(trim(adjustl(name12))== trim(adjustl(nameparent(ipar))))then
					JKP=ipar
					ILIG_TYPE=MET_TYPE(IPAR)! what type of metabolite is the template ligand?
					exit
				end if
			end do
			IF(JKP==0) go to 200

			ni=0
			do ia=1,naln
				i=mm2a(ia)
				if(i==0)go to 200
				ires=jcode(i)
				j=mm1a(ia)
				if(j==0)go to 200
				jres=icode(j,JKP)
				if(ires==jres)ni=ni+1
	    	end do
			ll=naln

! rotate template to target pocket
			do ia=1,naln
				i=mm2a(ia)
				j=mm1a(ia)
				do jj=1,3
	            	r_2(jj,ia)=ca(jj,i)
		      		r_1(jj,ia)=ca_par(jj,j,JKP)
				end do
			enddo
  	 		call u3b(w,r_1,r_2,ll,1,drms,u,t,ier) !u rotate r_1 to r_2
	     	armsd=dsqrt(drms/ll)
			if(armsd>3.5)go to 200

      		nh=nh_par(JKP)
       		do ih=1,nh
       			do j=1,3
                    x(j)=xhet(j,ih,JKP)
            	end do
            	xrot(1,ih)=t(1)+u(1,1)*x(1)+u(1,2)*x(2)+u(1,3)*x(3)
          		xrot(2,ih)=t(2)+u(2,1)*x(1)+u(2,2)*x(2)+u(2,3)*x(3)
            	xrot(3,ih)=t(3)+u(3,1)*x(1)+u(3,2)*x(2)+u(3,3)*x(3)
       		end do

			io=0
			itouch=0
			do i=1,mres
			do ih=1,nh
				dd=0.
				do j=1,3
					dd=dd+(xrot(j,ih)-ca(j,i))**2 ! refers to target protein
				enddo
!				if(dd.lt.25)then	! cannot have too many excluded volume overlaps
				if(dd.lt.9)then	! cannot have too many excluded volume overlaps
					io=io+1 
				end if
				if(dd.lt.30)then ! ligand must touch the target proten
					itouch=itouch+1
				end if		
			end do
			end do

			if(itouch <3) go to 200! rotated ligand must touch the protein
			if(io.gt.2*nh)then !cannot have too many ligand-protein overlaps
				go to 200
			end if
			IF(NI>6)NI=6
			
			NEX(NI,ILIG_TYPE)=NEX(NI,ILIG_TYPE)+1
			NN=NEX(NI,ILIG_TYPE)
			
			PV_MET(NI,NN,ILIG_TYPE)=pv
			
			NG=NG+1
			NLIG(ILIG_TYPE)=NLIG(ILIG_TYPE)+1
			TASS(ILIG_TYPE)=.true.
200		continue
close(591)
201	CONTINUE
	open(unit=59,Action='write',file='results_prec_0.05_'//trim(adjustl(targ))//'/'//trim(adjustl(namepdb(IK)))//'.summary')
!	write(6,*)'opening',namepdb(ik)//'.summary'
!	write(6,*)'NG=',NG
	NCS=0
	
	DO IM=1,NMET
		IF(TASS(IM))NCS=NCS+1
	END DO
	write(59,29)NCS
!	write(9,29)NCS !total number of ligands
	29     format(i4)
write(9,39)name12T,NCS
39 format(a12,1x,I4)
IF(NCS==0)CYCLE
NTOT_PREC=NTOT_PREC+1
NT=0
DO IM=1,NMET
	IF(TASS(IM))then
			precision=0.0
			ipsel=0
			NG=0
			DO NI=0,6
				IF(NEX(NI,IM)>0)then
					NG=NG+1
					DO NEX_I=1,NEX(NI,IM)
						ND=NEX(NI,IM)
						IF(ND>6)ND=6
						PV=PV_MET(NI,NEX_I,IM)
						IF(PV<=0.005)then	
							ip=1
						else
							ip=2
						end if
						if(prec(NI,ND,IP)>precision)then
							ipsel=ip
							NI_SEL=NI
							NEX_SEL=ND
							precision=prec(NI,ND,IP)
							rec_sel=rec(ni,ND,ip)
						end if
					end do
				end if
			ENDDO
			IF(NG>0)then
				NT=NT+1
			prec_met(NT)=precision
                        rec_met(NT)=rec_sel
                        anpv_met(NT)=anpv_sel
                        NI_D(NT)=NI_SEL
                        NEX_D(NT)=NEX_SEL
                        PV_D(NT)=PV_SEL(IPSEL)
                        MET_D(NT)=MET(IM)
				write(59,50)MET(IM),precision,rec_sel,NI_SEL,NEX_SEL,PV_SEL(IP)
				write(9,50)MET(IM),precision,rec_sel,NI_SEL,NEX_SEL,PV_SEL(IP)
				write(6,50)MET(IM),precision,rec_sel,NI_SEL,NEX_SEL,PV_SEL(IP)
				50 format(a3,1x,2(1f7.3,1x),2(i3,1x),1pd9.3)
			END IF
	END IF
	END DO
	write(9,*)'======='
!	write(6,*)namepdb(ik),NCS," no case
	write(6,*)'completed:',name12T
	close(59)
        do id=1,nt
                item(id)=id
	end do
 	do id=1,NT-1
                do jd=id+1,NT
                if(prec_met(id)<prec_met(jd))then
			ap=prec_met(id)
     			prec_met(id)=prec_met(jd)
                        prec_met(jd)=ap
                        it=item(id)
                        item(id)=item(jd)
                        item(jd)=it
                endif
                end do
        end do
        write(19,*)namePDB(IK),NT
        write(19,*)namePDB(IK),NT
        write(29,*)namePDB(IK),MIN(NT,5)
        write(6,*)namePDB(IK),MIN(NT,5)
        do iim=1,NT
                im=item(iim)
                if(iim<6)then
                        write(29,19)MET_D(im),prec_met(iim),rec_met(im),anpv_met(im),ni_D(im),nex_D(im),PV_D(im)
                        write(6,19)MET_D(im),prec_met(iim),rec_met(im),anpv_met(im),ni_D(im),nex_D(im),PV_D(im)
                end if
19		format(a3,1x,3(1f7.3,1x),2(i3,1x),1pd9.3)
        end do
        close(19)
        close(29)
2999 CONTINUE
close(9)

END
	


      subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
      integer ip(9), ip2312(4), i, j, k, l, m1, m, ier, n, mode
      PARAMETER(nmax=2000)
      double precision :: w(nmax), x(3, nmax), y(3, nmax), u(3, 3), t(3), rms, sigma
      double precision :: r(3, 3), xc(3), yc(3), wc, a(3, 3), b(3, 3), e0, e(3), e1, e2, e3, d, spur, det, cof, h, g, cth, sth, sqrth, p, tol, rr(6), rr1, rr2, rr3, rr4, rr5, rr6, ss(6), ss1, ss2, ss3, ss4, ss5, ss6, zero, one, two, three, sqrt3
      equivalence (rr(1), rr1), (rr(2), rr2), (rr(3), rr3), (rr(4), rr4),(rr(5), rr5), (rr(6), rr6), (ss(1), ss1), (ss(2), ss2), (ss(3), ss3), (ss(4), ss4), (ss(5), ss5), (ss(6), ss6), (e(1), e1), (e(2), e2), (e(3), e3)
      data sqrt3 / 1.73205080756888d+00 /
      data tol / 1.0d-2 /
      data zero / 0.0d+00 /
      data one / 1.0d+00 /
      data two / 2.0d+00 /
      data three / 3.0d+00 /
      data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
      data ip2312 / 2, 3, 1, 2 /
      wc = zero
	!!write(6,*)'n',n
      rms = 0.0
      e0 = zero
      do 1 i = 1, 3
      xc(i) = zero
      yc(i) = zero
      t(i) = 0.0
      do  110 j = 1, 3
      	d = zero
      	if (i .eq. j) d = one
      	u(i,j) = d
      	a(i,j) = d
     	r(i,j) = zero
	110 continue
1 	continue
      ier = -1
!**** DETERMINE CENTROIDS OF BOTH VECTOR SETS X AND Y
! 170 "rms.for"
      if (n .lt. 1) return 
! 172 "rms.for"
      ier = -2
      do 2 m = 1, n
      if (w(m) .lt. 0.0) return 
      wc = wc + w(m)
      do 210 i = 1, 3
      xc(i) = xc(i) + (w(m) * x(i,m))
      yc(i) = yc(i) + (w(m) * y(i,m))
	210 continue
	2	continue
      if (wc .le. zero) return 
      do 3 i = 1, 3
      xc(i) = xc(i) / wc
!**** DETERMINE CORRELATION MATRIX R BETWEEN VECTOR SETS Y AND X
! 182 "rms.for"
     yc(i) = yc(i) / wc
	 3 continue
! 184 "rms.for"
      do 4 m = 1, n
      do 41 i = 1, 3
      e0 = e0 + (w(m) * (((x(i,m) - xc(i)) ** 2) + ((y(i,m) - yc(i)) **2)))
! 187 "rms.for"
      d = w(m) * (y(i,m) - yc(i))
      do 42 j = 1, 3
! *** CALCULATE DETERMINANT OF R(I,J)
! 189 "rms.for"
    	r(i,j) = r(i,j) + (d * (x(j,m) - xc(j)))
		42 continue
		41 continue
		4 continue
! 191 "rms.for"
      det = ((r(1,1) * ((r(2,2) * r(3,3)) - (r(2,3) * r(3,2)))) - (r(1,2) * ((r(2,1) * r(3,3)) - (r(2,3) * r(3,1))))) + (r(1,3) * ((r(2,1)* r(3,2)) - (r(2,2) * r(3,1))))
! *** FORM UPPER TRIANGLE OF TRANSPOSED(R)*R
! 194 "rms.for"
      sigma = det
! 196 "rms.for"
      m = 0
      do 5 j = 1, 3
      do 510 i = 1, j
      m = m + 1
! **************** EIGENVALUES *****************************************
! *** FORM CHARACTERISTIC CUBIC  X**3-3*SPUR*X**2+3*COF*X-DET=0
! 200 "rms.for"
     	rr(m) = ((r(1,i) * r(1,j)) + (r(2,i) * r(2,j))) + (r(3,i) * r(3,j))
	510 continue
	5  continue
! 203 "rms.for"
      spur = ((rr1 + rr3) + rr6) / three
      cof = ((((((rr3 * rr6) - (rr5 * rr5)) + (rr1 * rr6)) - (rr4 * rr4)) + (rr1 * rr3)) - (rr2 * rr2)) / three
! 205 "rms.for"
      det = det * det
      do 6 i = 1, 3
    	e(i) = spur
	6	continue
!**** REDUCE CUBI! TO STANDARD FORM Y**3-3HY+2G=0 BY PUTTING X=Y+SPUR
! 208 "rms.for"
      if (spur .le. zero) goto 40
! 210 "rms.for"
      d = spur * spur
      h = d - cof
! *** SOLVE CUBIC. ROOTS ARE E1,E2,E3 IN DECREASING ORDER
! 212 "rms.for"
      g = (((spur * cof) - det) / two) - (spur * h)
! 214 "rms.for"
      if (h .le. zero) goto 8
      sqrth = dsqrt(h)
      d = ((h * h) * h) - (g * g)
      if (d .lt. zero) d = zero
      d = datan2(dsqrt(d),- g) / three
      cth = sqrth * dcos(d)
      sth = (sqrth * sqrt3) * dsin(d)
      e1 = (spur + cth) + cth
      e2 = (spur - cth) + sth
      e3 = (spur - cth) - sth
!.....HANDLE SPECIAL CASE OF 3 IDENTICAL ROOTS
! 224 "rms.for"
	  if(mode==0)then	
		go to 50
	  else 
		go to 10
	  end if
		
      !if (mode) 10, 50, 10
! *************** EIGENVECTORS *****************************************
! 226 "rms.for"
    8 if(mode==0)then	
		go to 50
	else
		go to 30
	end if
	!if (mode) 30, 50, 30
! 228 "rms.for"
   10 do 15 l = 1, 3, 2
      d = e(l)
      ss1 = ((d - rr3) * (d - rr6)) - (rr5 * rr5)
      ss2 = ((d - rr6) * rr2) + (rr4 * rr5)
      ss3 = ((d - rr1) * (d - rr6)) - (rr4 * rr4)
      ss4 = ((d - rr3) * rr4) + (rr2 * rr5)
      ss5 = ((d - rr1) * rr5) + (rr2 * rr4)
      ss6 = ((d - rr1) * (d - rr3)) - (rr2 * rr2)
      j = 1
      if (dabs(ss1) .ge. dabs(ss3)) goto 12
      j = 2
      if (dabs(ss3) .ge. dabs(ss6)) goto 13
   11 j = 3
      goto 13
   12 if (dabs(ss1) .lt. dabs(ss6)) goto 11
   13 d = zero
      j = 3 * (j - 1)
      do 14 i = 1, 3
      k = ip(i + j)
      a(i,l) = ss(k)
    d = d + (ss(k) * ss(k))
	14 continue
      if (d .gt. zero) d = one / dsqrt(d)
      do 151 i = 1, 3
    a(i,l) = a(i,l) * d
	151 continue
	15	continue
      d = ((a(1,1) * a(1,3)) + (a(2,1) * a(2,3))) + (a(3,1) * a(3,3))
      m1 = 3
      m = 1
      if ((e1 - e2) .gt. (e2 - e3)) goto 16
      m1 = 1
      m = 3
   16 p = zero
      do 17 i = 1, 3
      a(i,m1) = a(i,m1) - (d * a(i,m))
   p = p + (a(i,m1) ** 2)
   17	continue
      if (p .le. tol) goto 19
      p = one / dsqrt(p)
      do 18 i = 1, 3
    	a(i,m1) = a(i,m1) * p
	18 continue
      goto 21
   19 p = one
      do 20 i = 1, 3
      if (p .lt. dabs(a(i,m))) goto 20
      p = dabs(a(i,m))
      j = i
   20 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((a(k,m) ** 2) + (a(l,m) ** 2))
      if (p .le. tol) goto 40
      a(j,m1) = zero
      a(k,m1) = - (a(l,m) / p)
      a(l,m1) = a(k,m) / p
   21 a(1,2) = (a(2,3) * a(3,1)) - (a(2,1) * a(3,3))
      a(2,2) = (a(3,3) * a(1,1)) - (a(3,1) * a(1,3))
! ***************** ROTATION MATRIX ************************************
! 282 "rms.for"
      a(3,2) = (a(1,3) * a(2,1)) - (a(1,1) * a(2,3))
! 284 "rms.for"
   30 do 32 l = 1, 2
      d = zero
      	do 31 i = 1, 3
     		 b(i,l) = ((r(i,1) * a(1,l)) + (r(i,2) * a(2,l))) + (r(i,3) * a(3,l))
! 288 "rms.for"
   			d = d + (b(i,l) ** 2)
		31	continue
      if (d .gt. zero) d = one / dsqrt(d)
      do 321 i = 1, 3
    	b(i,l) = b(i,l) * d
		321 continue
	32 continue
      d = ((b(1,1) * b(1,2)) + (b(2,1) * b(2,2))) + (b(3,1) * b(3,2))
      p = zero
      do 33 i = 1, 3
      b(i,2) = b(i,2) - (d * b(i,1))
    p = p + (b(i,2) ** 2)
   33	continue
      if (p .le. tol) goto 35
      p = one / dsqrt(p)
      do 34 i = 1, 3
    b(i,2) = b(i,2) * p
	34 continue
      goto 37
   35 p = one
      do 36 i = 1, 3
      if (p .lt. dabs(b(i,1))) goto 36
      p = dabs(b(i,1))
      j = i
   36 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((b(k,1) ** 2) + (b(l,1) ** 2))
      if (p .le. tol) goto 40
      b(j,2) = zero
      b(k,2) = - (b(l,1) / p)
      b(l,2) = b(k,1) / p
   37 b(1,3) = (b(2,1) * b(3,2)) - (b(2,2) * b(3,1))
      b(2,3) = (b(3,1) * b(1,2)) - (b(3,2) * b(1,1))
      b(3,3) = (b(1,1) * b(2,2)) - (b(1,2) * b(2,1))
      do 39 i = 1, 3
      do 391 j = 1, 3
! ***************** TRANSLATION VECTOR *********************************
! 320 "rms.for"
    u(i,j) = ((b(i,1) * a(j,1)) + (b(i,2) * a(j,2))) + (b(i,3) * a(j,3))
	391 continue
	39 continue
   40 do 410 i = 1, 3
! ********************* RMS ERROR **************************************
! 323 "rms.for"
    t(i) = ((yc(i) - (u(i,1) * xc(1))) - (u(i,2) * xc(2))) - (u(i,3) * xc(3))
   410	continue
   50 do 51 i = 1, 3
      if (e(i) .lt. zero) e(i) = zero
   		e(i) = dsqrt(e(i))
   51 continue
      ier = 0
      if (e2 .le. (e1 * 1.0d-05)) ier = -1
      d = e3
      if (sigma .ge. 0.0) goto 52
      d = - d
      if ((e2 - e3) .le. (e1 * 1.0d-05)) ier = -1
   52 d = (d + e2) + e1
      rms = (e0 - d) - d
      if (rms .lt. 0.0) rms = 0.0
      return 
!.....END U3B...........................................................
!----------------------------------------------------------
!                       THE END
!----------------------------------------------------------
! 338 "rms.for"
      end

	
