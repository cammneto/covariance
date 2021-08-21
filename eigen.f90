!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  Programa para calcular a matriz covariância e  !!!!!!!
!!!!!!  seus respectivos autovalores a partir de um  !!!!!!!!
!!!!!!!!!  grupo de dados quaisquer  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  Autor:  Carlos Moreira de Melo Neto  !!!!!!!!!!!
!!!!!!!!!!!!!!!!  e-mail: cammneto@gmail.com !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  program eigen
  implicit none
   integer N, i
   real*8 covxx, covxy, covyy, xbar, ybar, x, y, x1, x2, x3, x4
   real L2, L1, D, C1, C2, C3, C4, C5, C6, C7, C8

   open(10,file='xy.dat')
   open(30,file='matcov.dat')
   open(40,file='médias.dat')
   open(50,file='autovalores.dat')
   open(60,file='autovetores.dat')
   open(70,file='confirma.dat')

     N=10 !Número de itens no grupo de dados
     xbar=0
     ybar=0

   do i=1,N
    read(10,FMT='(F4.1,2x,F4.1)') x, y
     xbar = xbar + x/N
     ybar = ybar + y/N
   end do

   rewind 10
   write (40,FMT='(A6,F5.2)'), 'xbar =', xbar
   write (40,FMT='(A6,F5.2)'), 'ybar =', ybar
   write (40,*), 'x-xbar y-ybar'
   do i=1,N
    read(10,FMT='(F4.1,2x,F4.1)') x, y
     covxx = covxx + ((x-xbar)**2)/(N-1)
     covyy = covyy + ((y-ybar)**2)/(N-1)
    write (40,FMT='(F5.2,2x,F5.2)'), x-xbar, y-ybar
   end do

   rewind 10
   do i=1,N
    read(10,FMT='(F4.1,2x,F4.1)') x, y
     covxy = covxy + ((x-xbar)*(y-ybar))/(N-1)
   end do

   write (30,*),'matriz covariância'
   write (30,*), covxx, covxy
   write (30,*), covxy, covyy

    D = sqrt((-covxx-covyy)**2 -4*(covxx*covyy-(covxy**2)))

    L1 = (covxx + covyy + D)/2

    L2 = (covxx + covyy - D)/2

   write (50,*),'Autovalores da matriz covariância'
   write (50,*), L1
   write (50,*), L2

!autovetores do autovalor L1
    x1 = covxy/sqrt((covxy**2)+(L1-covxx)**2)

    x2 = (L1-covxx)*x1/covxy

!autovetores do autovalor L2
    x3 =  covxy/sqrt((covxy**2)+(L2-covxx)**2)

    x4 = (L2-covxx)*x3/covxy

   write (60,*),'Autovetores para o autovalor =' , L1
   write (60,*), x1, x2
   write (60,*),'Autovetores para o autovalor =' , L2
   write (60,*), x3, x4

!confirmação para o autovalor L1
    C1 = covxx*x1 + covxy*x2
    C2 = covxy*x1 + covyy*x2
    C3 = x1*L1
    C4 = x2*L1
   write (70,*), 'Se as igualdades estão corretas encontramos autovalores e autovetores corretamente.'
   write (70,*), C1, '=', C3
   write (70,*), C2, '=', C4

!confirmação para o autovalor L2
    C5 = covxx*x3 + covxy*x4
    C6 = covxy*x3 + covyy*x4
    C7 = x3*L2
    C8 = x4*L2
   write (70,*), C5, '=',  C7
   write (70,*), C6, '=',  C8

end program
