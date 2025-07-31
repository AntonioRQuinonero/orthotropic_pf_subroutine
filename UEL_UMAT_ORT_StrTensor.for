      module kvisual
      implicit none
	  integer nelem, numsdv, nipt
	  parameter (nelem = 27500, numsdv=15, nipt=4)
	  !nelem = number of elements in an element layer
	  !numsdv = number of variables to be stored in sdvout (currently 15)
	  !ninpt = number of integration points in the element
	  
	  !sdvout(j,k,1) = Strain E11 in the PF element
	  !sdvout(j,k,2) = Strain E22 in the PF element
	  !sdvout(j,k,3) = Strain E12 in the PF element
	  !sdvout(j,k,4) = Strain E11 in the mechanical element
	  !sdvout(j,k,5) = Strain E22 in the mechanical element
	  !sdvout(j,k,6) = Strain E12 in the mechanical element
	  !sdvout(j,k,7) = PF variable phi1 in the PF element
	  !sdvout(j,k,8) = PF variable phi1 in the mechanical element
	  !sdvout(j,k,9) = PF variable phi2 in the PF element
	  !sdvout(j,k,10) = PF variable phi2 in the mechanical element
	  !sdvout(j,k,11) = empty
	  !sdvout(j,k,12) = Time instant
	  !sdvout(j,k,13) = Iteration number
	  !sdvout(j,k,14) = Fracture energy damage mechanism 1
	  !sdvout(j,k,15) = Fracture energy damage mechanism 2
	  

	  real*8 sdvout(nelem,nipt,numsdv)
	  
      save
      end module    


********************************************************************************
C UEL subroutine for a 2D plane stress full integration hex element
********************************************************************************
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)

	  use kvisual
      INCLUDE 'ABA_PARAM.INC'
	  	  
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
	 	  
  
	  parameter(gaussPoint=0.577350269189626d0,pi=3.14159265358979d0)

	  parameter(nintp=4, nodeDOF=6, ndim=2, ntens=3)
! 	  nintp = Number of integration points
!         nodeDOF = Number of degrees of freedom at a node
!	  ndim = Number of dimensions of the problem
!	  ntens = Size of the stress and strain vectors (plane stress=3)

	  dimension coordRef(4,2), dfunN(nnode,ndim), aInvMatJ(ndim,ndim), 
     1			Dmat(ntens,ntens),BmatPhi(ndim,nnode),
     2	        BmatPhiT(nnode,ndim), funN(nnode), dmat0(ntens,ntens),
     3			aMatJ(ndim,ndim), dN(nnode,ndim),force(ndim*nnode),
     4			stress(ntens),stran(ntens),gradphi1(ndim), gradphi2(ndim),
     6			dAMATRXPhi1(nnode,nnode),A3(nnode,nnode),
     7			damatrxphi1phi2(nnode,nnode),flux1(nnode), flux2(nnode),
     8			Cphi1(ntens,ntens), Cphi2(ntens,ntens),
     9			dAMATRXPhi2(nnode,nnode), Cphi1phi2(ntens,ntens),
     1		    StrTensor1(ndim,ndim), StrTensor2(ndim,ndim)


	  !Initialize sdvout
	  if (Time(2) .eq. 0.d0) then
		do i = 1,nelem
			do j = 1,nipt
				do k = 1,numsdv
					sdvout(i,j,k) = 0.d0
				enddo
			enddo
		enddo
	  endif
	  
	  !Reference element node coordinates
	  coordRef(1,1) = -1.
	  coordRef(1,2) = -1.
	  coordRef(2,1) = 1.
	  coordRef(2,2) = -1.
	  coordRef(3,1) = -1.
	  coordRef(3,2) = 1.
	  coordRef(4,1) = 1.
	  coordRef(4,2) = 1.

 	  !Material properties
	  E1 = props(1)
	  E2 = props(2)
	  xnu = props(3)
	  G12 = props(4)
	  xl1 = props(5)
	  Gc1 = props(6)
	  xl2 = props(7) 
	  Gc2 = props(8)
	  ang = props(9)*pi/180.d0

	  EnergiaEl = 0.d0 !Strain energy initialization
	  FractureEnergy1 = 0.d0 !Fracture energy 1 initialization
	  FractureEnergy2 = 0.d0 !Fracture energy 2 initialization
	  
          !Updating time and iteration variables

	  timez = sdvout(jelem,1,12)
	  if (timez .lt. time(2)) then
		sdvout(jelem,1,12) = time(2)
		sdvout(jelem,1,13) = 0
	  else
		sdvout(jelem,1,13) = sdvout(jelem,1,13) + 1
	  endif
	  stepiter = sdvout(jelem,1,13)
	 
	  !Element residuals and sitffness matrix initialization
	  do i = 1,nnode*nodeDof
		do j = 1,nnode*nodeDof
			AMATRX(i,j) = 0.d0
		enddo
		RHS(i,1) = 0.d0
	  enddo
	  
	  !Integration loop
	  do kintp=1,nintp
	  
!************** Shape functions and derivatives with respect local coordinates	  

		xi = coordRef(kintp,1)*gaussPoint
		eta = coordRef(kintp,2)*gaussPoint
		
		!Shape functions
		funN(3) = 0.25d0*(xi+1.)*(eta+1.)
		funN(4) = 0.25d0*(1.-xi)*(eta+1.)
		funN(2) = 0.25d0*(xi+1.)*(1.-eta)
		funN(1) = 0.25d0*(xi-1.)*(eta-1.)
		
		!Shape functions derivatives with respect to xi
		dfunN(3,1) = 0.25d0*(eta+1.)
		dfunN(4,1) = -0.25d0*(eta+1.)
		dfunN(2,1) = 0.25d0*(1.-eta)
		dfunN(1,1) = 0.25d0*(eta-1.)
		
		!Shape functions derivatives with respect to eta
		dfunN(3,2) = 0.25d0*(xi+1.)
		dfunN(4,2) = 0.25d0*(1.-xi)
		dfunN(2,2) = -0.25d0*(xi+1.)
		dfunN(1,2) = 0.25d0*(xi-1.)
		

!************** Jacobian (and inverse) matrix, determinant and shape function derivatives with respect global coordinates
		dxdxi = 0.d0
		dxdeta = 0.d0
		dydxi = 0.d0
		dydeta = 0.d0
	
		!Derivatives of the global coordiantes with respect local coordinates
		do k=1,NNODE		
			dxdxi = dxdxi + dfunN(k,1)*COORDS(1,k)
			dxdeta = dxdeta + dfunN(k,2)*COORDS(1,k)
			dydxi = dydxi + dfunN(k,1)*COORDS(2,k)
			dydeta = dydeta + dfunN(k,2)*COORDS(2,k)		
		enddo
	
		!Jacobian matrix and determinant
		aMatJ(1,1) = dxdxi
		aMatJ(1,2) = dydxi
		aMatJ(2,1) = dxdeta
		aMatJ(2,2) = dydeta
		
		detMatJ = aMatJ(1,1)*aMatJ(2,2)-aMatJ(1,2)*aMatJ(2,1)
		
		
		!Inverse Jacobian matrix
		aInvMatJ(1,1) = aMatJ(2,2) / detMatJ
		aInvMatJ(1,2) = - aMatJ(1,2) / detMatJ
		aInvMatJ(2,1) = - aMatJ(2,1) / detMatJ
		aInvMatJ(2,2) = aMatJ(1,1) / detMatJ
	  
!********************* B Matrix (PF problem)

		!Shape functions derivatives with respect global coordinates
		do k = 1, NNODE
			dN(k,1) = aInvMatJ(1,1)*dfunN(k,1) + aInvMatJ(1,2)*dfunN(k,2)
			dN(k,2) = aInvMatJ(2,1)*dfunN(k,1) + aInvMatJ(2,2)*dfunN(k,2)
		enddo
		
		!B matrix and transpose B
		do i = 1,ndim
			do j = 1,nnode
				BmatPhi(i,j) = dN(j,i)
				BmatPhiT(j,i) = dN(j,i)
			enddo
		enddo

!****************** Residual
	
	  !Se calcula phi y su gradiente a nivel punto de integración
	  	phi1 = 0.d0
		dphi1 = 0.d0
		gradphi1 = 0.d0
		phi2 = 0.d0
		dphi2 = 0.d0
		gradphi2 = 0.d0
		
	  do i = 1,nnode		
		phi1 = phi1 + funN(i)*u(1+nodeDof*(i-1))
		phi2 = phi2 + funN(i)*u(2+nodeDof*(i-1))
		
		dphi1 = dphi1 + funN(i)*du(1+nodeDof*(i-1),1)
		dphi2 = dphi2 + funN(i)*du(2+nodeDof*(i-1),1)
		do j = 1,ndim
			gradphi1(j) = gradphi1(j) + BmatPhi(j,i)*u(1+nodeDof*(i-1))
			gradphi2(j) = gradphi2(j) + BmatPhi(j,i)*u(2+nodeDof*(i-1))
		enddo
	  enddo
	  
	  !Actualizo phi mec tal y como hace molnar
	  if (stepiter .eq. 0.d0) then
		sdvout(jelem,kintp,7) = phi1-dphi1
		sdvout(jelem,kintp,9) = phi2-dphi2
	  else
		sdvout(jelem,kintp,7) = phi1
		sdvout(jelem,kintp,9) = phi2
	  endif
	  
	  !Se calculan funciones de degradación y geométricas. De momento, modelo AT2
	  !Atento a la forma de definir la degradación!!!!!
	  xk = 1.d-07 !Parámetro para rigidez residual en degradación completa
	  
	  ! DEGRADATION FUNCTION

	  g1 = 1.d0-phi1
	  g2 = 1.d0-phi2

	   
	  !GEOMETRIC FUNCTION
      alpha1 = phi1**2.d0
	  dalpha1 = 2.d0*phi1
	  ddalpha1 = 2.d0
	  cw1 = 2.d0
	  
	  alpha2 = phi2**2.d0
	  dalpha2 = 2.d0*phi2
	  ddalpha2 = 2.d0
	  cw2 = 2.d0
	  
	  !Matrices constitutivas y derivadas
	  do i=1,ntens
		do j=1,ntens
			Dmat(i,j) = 0.d0
			Dmat0(i,j) = 0.d0
		enddo
	  enddo
	  
	  b = 1.d0/(E1-E2*xnu**2.d0)
	  !Matriz constitutiva sin daño
	  Dmat0(1,1) = E1**2.d0*b
	  Dmat0(1,2) = E1*E2*xnu*b
	  Dmat0(2,2) = E1*E2*b
	  Dmat0(2,1) = Dmat0(1,2)
	  Dmat0(3,3) = G12
!	  
!	  !Matriz Constitutiva dañada
!	  do i = 1,NTENS
!		do j = 1,ntens
!			Dmat(i,j) = g*Dmat0(i,j)
!		enddo
!	  enddo
	    
	  !Lectura de deformaciones del problema mecánico
	  if (stepiter .eq. 0) then
		stran(1) = sdvout(jelem,kintp,4)
		stran(2) = sdvout(jelem,kintp,5)
		stran(3) = sdvout(jelem,kintp,6)
	  else
		stran(1) = sdvout(jelem,kintp,1)
		stran(2) = sdvout(jelem,kintp,2)
		stran(3) = sdvout(jelem,kintp,3)
	  endif

	  !Cáclulo de las derivadas de la matriz constituvia respecto variables PF
	  Cphi1 = 0.d0
	  Cphi1(1,1) = -2.d0*g1*Dmat0(1,1)
	  Cphi1(1,2) = -g2*Dmat0(1,2)
	  Cphi1(3,3) = -g2*Dmat0(3,3)
	  Cphi1(2,1) = Cphi1(1,2)
	  
	  Cphi2 = 0.d0
	  Cphi2(1,2) = -g1*Dmat0(1,2)
	  Cphi2(2,2) = -2.d0*g2*Dmat0(2,2)
	  Cphi2(3,3) = -g1*Dmat0(3,3)
	  Cphi2(2,1) = Cphi2(1,2)
	  
	  !Cálculo del tensor estructural en coordenadas globales
	  
	  !Tensor estructural 1
	  StrTensor1(1,1) = sin(ang)**2.d0
	  StrTensor1(1,2) = -cos(ang)*sin(ang)
	  StrTensor1(2,1) = -cos(ang)*sin(ang)
	  StrTensor1(2,2) = cos(ang)**2.d0
	  !Tensor estructural 2
	  StrTensor2(1,1) = cos(ang)**2.d0
	  StrTensor2(1,2) = cos(ang)*sin(ang)
	  StrTensor2(2,1) = cos(ang)*sin(ang)
	  StrTensor2(2,2) = sin(ang)**2.d0
	  
	  !Vectores base globales
	  
	  !Cálculo del residuo
	  !variable phi
	  
	  flux1 = 0.d0
	  flux2 = 0.d0
	
	  flux1 = (0.5d0*dot_product(stran,matmul(Cphi1,stran))*funN + 
     1			Gc1/(xl1*cw1)*(dalpha1*funN +
     2			2.d0*xl1**2.d0*matmul(BmatPhiT,matmul(StrTensor1,
     3			gradphi1))))*detMatJ
	 
	  flux2 = (0.5d0*dot_product(stran,matmul(Cphi2,stran))*funN + 
     1			Gc2/(xl2*cw2)*(dalpha2*funN +
     2			2.d0*xl2**2.d0*matmul(BmatPhiT,matmul(StrTensor2,
     3			gradphi2))))*detMatJ
	 	 
	  do i = 1,nnode
		RHS(1+nodeDof*(i-1),1) = RHS(1+nodeDof*(i-1),1) - flux1(i)
		RHS(2+nodeDof*(i-1),1) = RHS(2+nodeDof*(i-1),1) - flux2(i)
	  enddo
  	  
	  !Aprovecho para calcular la energía de fractura
	  gradMod1 = dot_product(gradphi1,matmul(StrTensor1,gradphi1))
	  gradMod2 = dot_product(gradphi2,matmul(StrTensor2,gradphi2))
	
	  FractureEnergy1 = FractureEnergy1 + (alpha1/cw1/xl1+xl1/cw1*
     1					gradMod1)*Gc1*detMatJ  

	  sdvout(jelem,kintp,14) = (alpha1/cw1/xl1+xl1/cw1*
     1					gradMod1)*Gc1*detMatJ  
	 
	  FractureEnergy2 = FractureEnergy2 + (alpha2/cw2/xl2+xl2/cw2*
     1					gradMod2)*Gc2*detMatJ

	  sdvout(jelem,kintp,15) = (alpha2/cw2/xl2+xl2/cw2*
     1					gradMod2)*Gc2*detMatJ
	 
!*************** Element stiffness matrix

	  Cphi1phi2 = 0.d0
	  Cphi1phi2(1,2) = Dmat0(1,2)
	  Cphi1phi2(2,1) = Dmat0(1,2)
	  Cphi1phi2(3,3) = Dmat0(3,3)
	  
	  !Matriz auxiliar para todos los cálculos posteriores
	  do i = 1,nnode
		do j = 1,nnode
			A3(i,j) = funN(i)*funN(j)
		enddo
	  enddo
	  
	  !Kphi1phi1
	  dAmatrxPhi1 = (Dmat0(1,1)*stran(1)**2.d0*A3 + 
     1				Gc1/(cw1*xl1)*(ddalpha1*A3 + 
     2				2.d0*xl1**2.d0*matmul(BmatPhiT,matmul(StrTensor1,
     3				BmatPhi))))*detMatJ
	  
	  !Kphi2phi2
	  dAmatrxPhi2 = (Dmat0(2,2)*stran(2)**2.d0*A3 + 
     1				Gc2/(cw2*xl2)*(ddalpha2*A3 + 
     2				2.d0*xl2**2.d0*matmul(BmatPhiT,matmul(StrTensor2,
     3				BmatPhi))))*detMatJ
	 
	  !Kphi1phi2
	  dAmatrxPhi1Phi2 = (0.5d0*dot_product(stran,matmul(Cphi1phi2,stran))*
     1					A3)*detMatJ
	  
	  
	  do i = 1,nnode
		do j = 1,nnode
		
		AMATRX(1+nodeDof*(i-1),1+nodeDof*(j-1)) = 
     1					AMATRX(1+nodeDof*(i-1),1+nodeDof*(j-1)) + 
     2					dAmatrxPhi1(i,j)
	 
		AMATRX(1+nodeDof*(i-1),2+nodeDof*(j-1)) = 
     1					AMATRX(1+nodeDof*(i-1),2+nodeDof*(j-1)) + 
     2					dAmatrxPhi1Phi2(i,j)
		
		AMATRX(2+nodeDof*(i-1),1+nodeDof*(j-1)) = 
     1					AMATRX(2+nodeDof*(i-1),1+nodeDof*(j-1)) + 
     2							dAmatrxPhi1Phi2(i,j)
	 
		AMATRX(2+nodeDof*(i-1),2+nodeDof*(j-1)) = 
     1					AMATRX(2+nodeDof*(i-1),2+nodeDof*(j-1)) + 
     2					dAmatrxPhi2(i,j)
		
		enddo
	  enddo

!************** Actualización de variables externas sdvout salida PF

	  sdvout(jelem,kintp,1) = stran(1)
	  sdvout(jelem,kintp,2) = stran(2)
	  sdvout(jelem,kintp,3) = stran(3)
		
	  
	  enddo
	  
 	  ResEnergy = 0.d0 !Energía "absorbida" por las restricciones
	  
	  !IRREVERSIBILIDAD
	   !Calculamos los residuos y Tanget stiffness matrix cuando salimos de las restricciones
	  do k = 1,nnode
	  
		!Comprobamos condiciones en phi1 en cada nodo
		if (du(1+nodeDof*(k-1),1) .lt. 0.d0 .or. 
     1	u(3+nodeDof*(k-1)) .gt. 0.d0) then
			RHS(1+nodeDof*(k-1),1) = RHS(1+nodeDof*(k-1),1) +
     1								 u(3+nodeDof*(k-1))
			RHS(3+nodeDof*(k-1),1) = RHS(3+nodeDof*(k-1),1) + 
     1								 du(1+nodeDof*(k-1),1)
			AMATRX(1+nodeDof*(k-1),3+nodeDof*(k-1)) = -1.d0
			AMATRX(3+nodeDof*(k-1),1+nodeDof*(k-1)) = -1.d0
		else
			AMATRX(3+nodeDof*(k-1),3+nodeDof*(k-1)) = 1.d0
		endif
		
		!Comprobamos condiciones en phi2 en cada nodo
		if (du(2+nodeDof*(k-1),1) .lt. 0.d0 .or. 
     1	u(4+nodeDof*(k-1)) .gt. 0.d0) then
			RHS(2+nodeDof*(k-1),1) = RHS(2+nodeDof*(k-1),1) + 
     1								 u(4+nodeDof*(k-1))
			RHS(4+nodeDof*(k-1),1) = RHS(4+nodeDof*(k-1),1) +
     1								 du(2+nodeDof*(k-1),1)
			AMATRX(2+nodeDof*(k-1),4+nodeDof*(k-1)) = -1.d0
			AMATRX(4+nodeDof*(k-1),2+nodeDof*(k-1)) = -1.d0
		else
			AMATRX(4+nodeDof*(k-1),4+nodeDof*(k-1)) = 1.d0
		endif
	  enddo
	  
	!LIMITACIÓN PHI<1
		   !Calculamos los residuos y Tanget stiffness matrix cuando salimos de las restricciones
	  do k = 1,nnode
	  
		!Comprobamos condiciones en phi1 en cada nodo
		if (u(1+nodeDof*(k-1)) .gt. 0.9999d0 .or. 
     1	u(5+nodeDof*(k-1)) .gt. 0.d0) then
			RHS(1+nodeDof*(k-1),1) = RHS(1+nodeDof*(k-1),1) -
     1								 u(5+nodeDof*(k-1))
			RHS(5+nodeDof*(k-1),1) = 1.d0 - u(1+nodeDof*(k-1))
			AMATRX(1+nodeDof*(k-1),5+nodeDof*(k-1)) = 1.d0
			AMATRX(5+nodeDof*(k-1),1+nodeDof*(k-1)) = 1.d0
		else
			AMATRX(5+nodeDof*(k-1),5+nodeDof*(k-1)) = 1.d0
		endif
		
		!Comprobamos condiciones en phi2 en cada nodo
		if (u(2+nodeDof*(k-1)) .gt. 0.9999d0 .or. 
     1	u(6+nodeDof*(k-1)) .gt. 0.d0) then
			RHS(2+nodeDof*(k-1),1) = RHS(2+nodeDof*(k-1),1) - 
     1								 u(6+nodeDof*(k-1))
			RHS(6+nodeDof*(k-1),1) = 1.d0 - u(2+nodeDof*(k-1))
			AMATRX(2+nodeDof*(k-1),6+nodeDof*(k-1)) = 1.d0
			AMATRX(6+nodeDof*(k-1),2+nodeDof*(k-1)) = 1.d0
		else
			AMATRX(6+nodeDof*(k-1),6+nodeDof*(k-1)) = 1.d0
		endif
	  enddo
	 
	  
	  !Guardo la energía de fractura debido a phi como Creep dissipation
	  ENERGY(3) = FractureEnergy1 + FractureEnergy2
	 
	  ! print*, 'AMATRX'
	  ! do i=1,24
	  ! print("(24F6.2)"), AMATRx(i,1:24) 
	  ! enddo
      RETURN
      END	  

  
	      
************************************************************************
C 		UMAT sub
************************************************************************

	  SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
	  use kvisual !Módulo que contiene sigout para conectar UEL con UMAT
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

	  dimension c0(ntens,ntens)
	  
      integer kelem
	  
	  !Inicializamos DDSDDE
	  do i = 1, ntens
		do j = 1,ntens
			ddsdde(i,j) = 0.d0
			c0(i,j) = 0.d0
		enddo
	  enddo

	  !Lectura propiedades mecánicas
	  E1 = props(1)
	  E2 = props(2)
	  xnu = props(3)
	  G12 = props(4)
	  
      kelem=int(noel-nelem) !Se asume siempre que los elementos a utilizar en UMAT
	  !van numerados después de los elementos para la UEL
	 

	  !Lectura de los valores de phi
	  stepiter = sdvout(kelem,1,13)
	  if (stepiter .eq. 0) then
		phi1 = sdvout(kelem,npt,7)
		phi2 = sdvout(kelem,npt,9)
	  else
		phi1 = sdvout(kelem,npt,8)
		phi2 = sdvout(kelem,npt,10)
	  endif

	  g1 = max(1.d0-phi1,1d-7)
	  g2 = max(1.d0-phi2,1d-7)
	  ! g1 = 1.d0-phi1
	  ! g2 = 1.d0-phi2
	  
	  !Cálculo DDSDDE plane stress
	  a = 1.d0/(E1-E2*xnu**2.d0)

	  c0(1,1) = E1**2.d0*a
	  c0(1,2) = E1*E2*xnu*a
	  c0(2,2) = E1*E2*a
	  c0(2,1) = c0(1,2)
	  c0(3,3) = G12
	  	  	  
	  ddsdde(1,1) = g1**2.d0*c0(1,1)
	  ddsdde(1,2) = g1*g2*c0(1,2)
	  ddsdde(2,2) = g2**2.d0*c0(2,2)
	  ddsdde(3,3) = g1*g2*c0(3,3)
	  ddsdde(2,1) = ddsdde(1,2)


	  !Actualizo stresses
	  stress = matmul(ddsdde,stran+dstran)
		
	  !Energía de deformación degradada
	  sse = 0.5d0*dot_product(stress,stran+dstran)
	  
	  !Variable energía deformación
	  ! H = 0.5d0*dot_product(matmul(c0,stran+dstran),stran+dstran)
	  
	  !Irreversibilidad
	  ! H = max(H,sdvout(kelem,npt,10))
	  
	  !Actualizo las variables de comunicación

	  sdvout(kelem,npt,4) = stran(1) + dstran(1)
	  sdvout(kelem,npt,5) = stran(2) + dstran(2)
	  sdvout(kelem,npt,6) = stran(3) + dstran(3)
	  sdvout(kelem,npt,8) = phi1
	  sdvout(kelem,npt,10) = phi2
	  

      do k1 = 1, numsdv-3
        statev(k1) = sdvout(kelem,npt,k1)
      end do 
		
		statev(11) = sdvout(kelem,npt,14)
		statev(12) = sdvout(kelem,npt,15)
		
      RETURN
      END
