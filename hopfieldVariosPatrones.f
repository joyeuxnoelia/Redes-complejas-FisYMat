      program hopfieldVariosPatrones
      
      implicit none
    

      integer i,j,k,L,contador(105),mc(105),q,total_mc,neuronasDef,bucle
      integer nPatrones, iteracion
      parameter (nPatrones=3)
      parameter(L=30)
      real*8 T,valor,epsilon,dran_u,p,aleat,def
      real*8 a(105,nPatrones),peso(105,L,L,L,L)
      real*8 umbral(105,L,L),varH(105)
      real*8 solapamiento(101,nPatrones),spromedio(101,nPatrones)
      integer*4 i_dran,n,m,patron(L,L,1:nPatrones)
      integer*4 s(105,L,L),sini(L,L)
      logical deformacion,estudiarSolapamientoPromedio

      open(14,file='patron1.dat')
      open(15,file='patron2.dat')
      open(16,file='patron3.dat')
      
      open(10,file='patronInicial.dat') !guarda solo el valor inicial
      open(11,file='avancePatron.dat') !guarda la variacion de la matriz
      open(2,file='patronFinal.dat')
      
      open(3,file='solapamiento1.dat')
      open(4,file='solapamiento2.dat')
      open(5,file='solapamiento3.dat')

      open(21,file='solapPromedio_P1_def0_3.dat')
      open(22,file='solapPromedio_P2_def0_3.dat')
      open(23,file='solapPromedio_P3_def0_3.dat')
      
      call dran_ini(18635112)
      
!************************************************************************
!     temperatura
      T=0.0001d0

      total_mc=100
      
      estudiarSolapamientoPromedio=.false.
      
      deformacion=.false.
      def=0.3d0
      m=2       !patron a deformar
!**************************************************************************
      
      !     Grabamos el patron

      do k=1,nPatrones
         do i=1,(L)
            read(13+k,*) (patron(i,j,k), j=1,L)
         enddo
      enddo

      !     Generación aleatoria de la matriz de 0 y 1 en caso de partir de una situación aleat.

         if(deformacion.eqv..false.)then
            
            do i=1,L
               do j=1,L
                  aleat=i_dran(2)
                  if(aleat.eq.1)then
                     sini(i,j)=1
                  else
                     sini(i,j)=0
                  endif
               enddo
            enddo
            
         else
            
!     Generacion de la matriz partiendo del patron deformado. 

           
            
            do i=1,L
               do j=1,L
                  sini(i,j)=patron(i,j,m)
               enddo
            enddo
            
            
            neuronasDef=int(def*L*L) 
            
         
            do k=1, neuronasDef
               i=i_dran(L)
               j=i_dran(L)
               
               sini(i,j)=1-sini(i,j)
            enddo
         
            
         endif
        
!     Se guarda en un archivo la matriz inicial obtenida al reves
         do i=1,L
            write(10,*) (s(q,L+1-i,j), j=1,L)
            write(11,*) (s(q,L+1-i,j), j=1,L)
         enddo
               
         write(11,*)
         write(11,*)
      


      if(estudiarSolapamientoPromedio.eqv..true.) then
         iteracion=100
      else
         iteracion=1
      endif
            
      do q=1,iteracion          ! para diferentes T
         
         if(estudiarSolapamientoPromedio.eqv..true.) then
            T=0.0001*10000.d0**(dfloat(q)/100.d0)
         else
            T=T
         endif
        
         contador(q)=0
        
         spromedio=0.d0
         solapamiento=0.d0

!     pasos de montecarlo

        
         mc(q)=1
         
!     valor inicial matriz de espines
         do i=1,L
            do j=1,L
               s(q,i,j)=sini(i,j)
            enddo
         enddo


!     Para calcular los pesos sinapticos, calculamos el promedio de neuronas activadas
         a=0.d0

         do k=1,nPatrones
            do i=1, L
               do j=1,L
                  a(q,k)=a(q,k)+1.d0*patron(i,j,k)
               enddo
            enddo
         
            a(q,k)=a(q,k)/(1.d0*(L*L))
 
         enddo

!     pesos sinapticos

         peso=0.d0
         
         do i=1,L
            do j=1,L
               do m=1,L
                  do n=1,L
                     if((i.eq.m).AND.(j.eq.n)) then
                        peso(q,i,j,m,n)=0.d0
                       
                     else
                        do k=1, nPatrones
                           peso(q,i,j,m,n)=peso(q,i,j,m,n)+(1.d0/
     $                          (1.d0*L*L))*(1.d0*patron
     $                          (i,j,k)-a(q,k))*(1.d0*patron(m,n,k)-
     $                          a(q,k))
                        enddo
                        
                     endif

                  enddo
               enddo
            enddo
         enddo
         
         
!     umbral 
         umbral=0.d0

         do i=1,L
            do j=1,L
               do m=1,L
                  do n=1,L
                     umbral(q,i,j)=umbral(q,i,j)+0.5d0*peso(q,i,j,m,n)
                     
                  enddo
               enddo
            enddo
         enddo
         
         
         
         do bucle=1,L*L*total_mc
!     se toma un punto al azar de la red
            
            n=i_dran(L)
            m=i_dran(L)


!     calcular h
            varH(q)=0.d0

            do i=1,L
               do j=1,L
                  varH(q)=varH(q)+peso(q,n,m,i,j)*s(q,i,j)
               enddo
            enddo

            varH(q)=(-varH(q)+umbral(q,n,m))*(1.d0-2.d0*s(q,n,m))
            
            
!     se evalúa p
            
            valor=exp(-(varH(q)/T))
            
            if(valor.lt.1) then
               p=valor
            else
               p=1.d0
            endif
            
            
!     cambio del spin (o no)
            
            epsilon=dran_u()
            
            
            if(epsilon.lt.p)then
               s(q,n,m)=1-s(q,n,m)
            endif
            
!     copiar los valores en un doc para cada paso
            contador(q)=contador(q)+1
            
            if(contador(q).eq.(L**2))then
       
               do i=1,L
                  write(11,*) (s(q,L+1-i,j), j=1,L)
               enddo
               
               write(11,*)
               write(11,*)

!     calculamos el solapamiento para cada paso
               solapamiento=0.d0
              

               do m=1,nPatrones
                  do i=1,L
                     do j=1,L
                        solapamiento(q,m)=solapamiento(q,m)+(1.d0*
     $                       patron(i,j,m)-a(q,m))*(1.d0*s(q,i,j)-
     $                       a(q,m))
                     enddo
                  enddo
       
                  solapamiento(q,m)=solapamiento(q,m)/(a(q,m)*(1.d0-
     $                 a(q,m))*1.d0*L*L)
                  write(2+m,*)mc(q), solapamiento(q,m)

!     solapamiento promedio
                  if(mc(q).ge.15) then
                     spromedio(q,m)=spromedio(q,m)+solapamiento(q,m)
                  endif
                  
               enddo

               mc(q)=mc(q)+1
              
               contador(q)=0
            endif
            
         enddo                  !bucle 

!     copio el resultado final a cada T en un archivo a parte

         do i=1,L
            write(2,*) (s(q,L+1-i,j), j=1,L)
         enddo

         write(2,*)
         write(2,*)

         do m=1,npatrones
            write(20+m,*)T,spromedio/86.d0
         enddo
         
      enddo                     !diferentes T

      
      close(14)
      close(15)
      close(16)
      
      close(10)
      close(11)
      close(2)
      
      close(3)
      close(4)
      close(5)

      close(21)
      close(22)
      close(23)
      

      stop
      end program
      
      include 'dranxor2_new.f'
