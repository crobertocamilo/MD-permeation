
!CITATION: If you use MD-permeation in your research, we ask that you cite the following article:
!C.R.S. Camilo, J.R. Ruggiero, and A.S. de Araujo. 2021. A method for detection of permeation events in Molecular Dynamics simulations of lipid bilayers. bioRxiv doi:10.1101/2021.01.20.427278

      program MD_permeation

      implicit none
     
      real xlp,ylp,zlp,vxlp,vylp,vzlp
      real x,y,z,vx,vy,vz
      real lsupp,linfp
      real top,lower,thick,center,t,ztop,zlower
      real t0,tn,ti
      real z0,z1,d1,d2,Ld,sim_time
      real tiP(5000),tfP(5000),dtP      !5000 are the maximum number of permeation events in the analyzed trajectory
      real xp,yp,zp,vxp,vyp,vzp

      character*7 molecule,atom
      character*80 cab
      character*4 residuo
      character*1 resp
                                        
      integer ps,nframes,t_cut,residP(5000),count_ev(50000)
      integer nro,natom,c1,c2,c3,c4,c5,c6,c7,c8,c9
      integer nmol,nwater,ntop,nlower,delta
      integer*4 k,l,m,n,gap_res,k0,cont,nevent,tf
      integer*4 contP(5000)            
                                      
      real Lx,Ly,Lz(50000)

      integer layer(10000,50000)        !layer's dimensions must be: layer(number of water molecules, number of frames)
                           

c     ******************************************************************
c                         INPUT FILES
c     ******************************************************************
      open(10,file='OW.gro',status='old')
      open(20,file='P.gro',status='old')
c     ******************************************************************
c                         OUTPUT FILES
c     ******************************************************************
      open(11,file='Permeation_events.txt',status='unknown')
      open(12,file='Thickness.txt',status='unknown')
c      open(13,file='Type2-jumps.txt',status='unknown')
      open(14,file='Checkpoint_delta.txt',status='unknown')
c     ******************************************************************

      write(*,*)"PERMEATION EVENTS ANALYSIS"
      write(*,*)
      write(*,*)"---Set the input parameters---"
      write(*,*)
      write(*,*)"Estimative for bilayer center (nm):"
      read(*,*)center

      ti=0   !If the initial time is not equal to zero, the 'Checkpoint_delta.txt' file from the previous analyses is required!

      write(*,*)"Simulation time (ns):"
      read(*,*)sim_time
      write(*,*)"Trajectory precision (ps):"
      read(*,*)ps

      nframes=sim_time*1000/ps+1
      write(*,*)
      write(*,*)'Number of trajectory frames to analyze:',nframes
      write(*,*)

      write(*,*)"Number of the first water residue:"
      read(*,*)gap_res
      write(*,*)

      gap_res=gap_res-1

      t_cut=0.0  !Treatment for type2+ jumps. Only accept events with a duration equal or greater than t_cut (in ps). 
	         !For a 20 ps trajectory precision, set t_cut=0.0 (all events will be acepted).
	         !For a 100ps trajectory precision, we recommend t_cut=0.4 or greater.



      write(12,*)"t"," ztop"," zlower"," thickness",
     &        " center"," ztop-center"," zlower-center" !'Thickness.txt' header: column names 
      write(12,*)"t(ns)"," ztop(nm)"," zlower(nm)"," thickness(nm)",
     &        " center(nm)"," top(A)"," lower(A)"       !'Thickness.txt' header: column comments
 
      write(*,*)'Reading the trajectories...'     

      do k=1,nframes 

       t=k
       t=ti+(t-1)/(1000/ps)

       ntop=0
       nlower=0
       ztop=0.0
       zlower=0.0

       read(20,*)cab

       read(20,*)natom

       do l=1,natom

          read(20,*)molecule,atom,nro,xlp,ylp,zlp

          if (zlp.GT.center) then

            ntop=ntop+1
            ztop=ztop+zlp

            else

            nlower=nlower+1
            zlower=zlower+zlp

          endif
       enddo   
       read(20,*)Lx,Ly,Ly

       top=ztop/ntop
       lower=zlower/nlower
      
       thick=top-lower
       center=lower+thick/2

       write(12,*)t,top,lower,thick,center,
     &             10*(top-center),10*(lower-center)

       lsupp=center+0.25*thick
       linfp=center-0.25*thick
      
   
       read(10,*)cab

       read(10,*)natom

       do l=1,natom

          z0=z    

          read(10,*)molecule,atom,nro,x,y,z,
     &              vx,vy,vz

c         Layer I: Above (outside) the bilayer 
          if (z.GT.top) then
             layer(l,k)=1
          endif

c         Layer V: Below (outside) the bilayer 
          if (z.LT.lower) then
	     layer(l,k)=5
          endif

c         Layer II: 0-25% inside the bilayer (top leaflet)
          if ((z.LE.top).and.(z.GT.lsupp)) then
             layer(l,k)=2
          endif       
      
c         Layer IV: 0-25% inside the bilayer (lower leaflet)
          if ((z.LT.linfp).and.(z.GE.lower)) then
             layer(l,k)=4
          endif   

c         Layer III: Center of the bilayer (hydrophobic nucleus)
          if ((z.LE.lsupp).and.(z.GE.linfp)) then
             layer(l,k)=3
          endif   

       enddo


c       write(*,*)"Reading: t =",t
     
       read(10,*)Lx,Ly,Lz(k) 
         
      enddo


      write(*,*)'Analyzing the data...'

      do l=1,natom
 
c          write(*,*)"Analysing: Molecule",l

c          delta=0
          cont=0
 
          if(ti.eq.0) then   !If the simulation is not a continuation (ti=0), the program will disregards the initial 5 ns of simulation. 
                k0=(5000/ps)-1
                delta=0
             else
                k0=0
                read(14,*)molecule,delta               
          endif 
          
          do k=(k0+1),nframes  

             if(abs(layer(l,k)-layer(l,k-1)).eq.1) then   
                delta=delta+(layer(l,k)-layer(l,k-1))
             endif


	     !Type2-jumps treatment.	

c             if(abs(layer(l,k)-layer(l,k-1)).eq.2) then   !Saltos do tipo 2 só são aceitos se a dist. direta (d1) for menor do que a dist. passando pela imagem (d2).
c
cc                z0=z(l,k-1)
c                z1=z
c                Ld=(Lz(k)+Lz(k-1))/2
c
c                d1=abs(z1-z0)  !Distância direta entre a posição da partícula em t e t + detal_t
c
c                if(z1.lt.z0)then !Distância considerando que a partícula saiu para a imagem - CPC
c                      d2=(Ld-z0)+z1
c                   else
c                      d2=z0+(Ld-z1)
c                endif
c
c                if(d1.le.d2) then
c                   delta=delta+(layer(l,k)-layer(l,k-1))
c                endif
c             endif

          
             if(abs(delta).eq.4) then

                t0=k0
                tn=k
                t0=ti+(t0-1)/(1000/ps)
                tn=ti+(tn-1)/(1000/ps)
                
                if((tn-t0).ge.t_cut)then

                   if(nevent.eq.0)then      !'Permeation_events.txt' header
                      write(11,*)"Permeation events of water molecules."
                      write(11,*)'Resid',' Initial_t(ns)',
     &                           ' Final_t(ns)',' dt(ns)'
                   endif
           
                   write(11,*)l+gap_res,t0,tn,tn-t0  

                   nevent=nevent+1

                   
                endif
             endif


             if(abs(delta).eq.3) then
                      
                if(((layer(l,k-1).eq.2).and.(layer(l,k).eq.5)).or.
     &             ((layer(l,k-1).eq.4).and.(layer(l,k).eq.1)).or.
     &             ((layer(l,k0).eq.5).and.(layer(l,k0+1).eq.2)).or.
     &             ((layer(l,k0).eq.1).and.(layer(l,k0+1).eq.4))) then

                t0=k0
                tn=k
                t0=ti+(t0-1)/(1000/ps)
                tn=ti+(tn-1)/(1000/ps)

                if((tn-t0).ge.t_cut)then

                   if(nevent.eq.0)then
                      write(11,*)"Permeation events of water molecules."
                      write(11,*)'Resid',' Initial_t(ns)',
     &                           ' Final_t(ns)',' dt(ns)'
                   endif
                                  
                   write(11,*)l+gap_res,t0,tn,(tn-t0)

                   nevent=nevent+1
                   
                endif

                endif
             endif

             if((layer(l,k).eq.1).or.(layer(l,k).eq.5)) then
                
                delta=0
                k0=k

             endif
     
                  !Recording all type2-jumps along the simulation
c             if(abs(layer(l,k)-layer(l,k-1)).eq.3) then
c                 write(13,*)l+gap_res,ti+(k-1)/(1000/ps)  
c             endif
          
          enddo

          close(30)
  
          write(residuo,90)l+gap_res
          write(14,*)residuo//'SOL',delta

      enddo

      write(*,*)
      write(*,*)nevent,' permeation events found.'
      write(*,*)
      write(*,*)

90    format(I4)      
          
      close(10)
      close(11)
      close(12)
c      close(13)
      close(14)
      close(20)

      
      open(11,file='Permeation_events.txt',status='old')
      open(16,file='Counting_events.txt',status='unknown')

      if(nevent.ne.0)then
	    write(16,*)"Permeation events over simulation time."
            write(16,*)'t(ns)',' Number of events'
         else
	    go to 404
      endif

      do l=1,nframes
         count_ev(l)=0
      enddo

      read(11,*)cab    !Reading header from "Permeation_events.txt" 
      read(11,*)cab

      do l=1,nevent

         read(11,*)residP(l),tiP(l),tfP(l),dtP

         do m=1,nframes

            t=m
            t=ti+(t-1)/(1000/ps)
        
            if(t.ge.tiP(l))then
               count_ev(m)=count_ev(m)+1
            endif
         enddo
      enddo

      do m=1,nframes
         t=m
         t=ti+(t-1)/(1000/ps)
         write(16,*)t,count_ev(m)
      enddo   


      write(*,*)"Do you want to map each event's individual trajectory?"
      write(*,*)"This may take several minutes! (y/n)?"
      read(*,*)resp
      write(*,*)


      if((resp=="y").or.(resp=="yes"))then

         do l=1,(natom+gap_res)
            contP(l)=0      
         enddo


         do l=1,nevent

            write(residuo,90)residP(l)   

            write(*,*)"Mapping...",l," from ",nevent,' (Residue: ',
     &         residuo,')'
                   
            if(contP(residP(l)).eq.0) then 

               contP(residP(l))=1   
               open(30,file=residuo//'SOL.txt',status='unknown')
               write(30,*)"t(ns)"," x"," y"," z"," vx"," vy"," vz"  !"residuoSOL.tx" header: column names
               write(30,*)residuo//"SOL ",residuo//"SOL ",          !"residuoSOL.tx" header: column comments
     &            residuo//"SOL ",residuo//"SOL ",residuo//"SOL ",
     &            residuo//"SOL ",residuo//"SOL"

            endif


            open(10,file='OW.gro',status='old')
            open(12,file='Thickness.txt',status='old')
 
            read(12,*)cab    !Reading header from "Thickness.txt"
            read(12,*)cab

            tf=(tfP(l)-ti)*(1000/ps)+1  

            do k=1,tf 

               read(12,*)t,ztop,zlower,thick,center  

               read(10,*)cab
               read(10,*)natom
               
               do n=1,natom

                  read(10,*)molecule,atom,nro,xp,yp,zp,vxp,vyp,vzp

                  if((n+gap_res).eq.residP(l))then
                     if((t.ge.tiP(l)).and.(t.le.tfP(l)))then   

                         write(30,*)t,10*xp,10*yp,            !The individual trajectories are centralized in relation to the center of the bilayer.
     &                          10*(zp-center),vxp,vyp,vzp
                     endif
                  endif
              
               enddo

c               write(*,*)"Event:",l,"Reading: t =",t
     
               read(10,*)Lx,Ly,Lx

            enddo

            close(10)
            close(12)

            if(residP(l).ne.residP(l+1))then  
               close(30)

            else
               write(*,*)'This residue has more than one event.'
            endif
 
         enddo

c	 write(*,*)"Plot 'Thickness.txt' and 'SOL.txt' file(s) at the same graph to visualize."
c    	 write(*,*)

      endif


404   close(10)
      close(20)
      close(11)
      close(16)


      write(*,*)'Done.'
      read(*,*)

      stop
      end


