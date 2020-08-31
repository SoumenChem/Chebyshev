C
      subroutine prop_chebyshev_poly(params,tt,it,dt,scalefock,tdfock,
     $     nmats, g_zfock, energies, g_zdens_mo,g_zdens_ao)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      integer, intent(in)           :: it
      double precision, intent(in)  :: tt, dt      !current time
      logical, external             :: tdfock  !external function for computing F(t)
      integer, intent(in)           :: nmats   !num of matrices (1 for CS and SO, 2 for OS)


C     == In/out ==
      integer, intent(in) :: g_zdens_mo(nmats)        ! in = P'(t); out = P'(t+dt)
      integer, intent(in) :: g_zfock(nmats)           ! F'(t) (MO basis); update at end

      
C     == Outputs ==
      type(rt_energies_t) energies


C     == Parameters ==
      character(*), parameter     :: pname =
     $     "prop_chebysev_poly: "
      integer, parameter          :: ninter_max = 200
      

C     == Variables ==
      integer me
      logical ok
      integer nmats2, nmats3
      integer imat
      integer it
      character amat(nmats)

      logical converged
      integer iinter
      double precision elapsed
      double precision scalefock(nmats)
      double precision ttdt     ! t+dt
      character*256 outstring
      integer num_same
      integer levals(nmats),ievals(nmats)

C     (simple way: have two GA handles and allocate one or two arrays, depending)
      integer g_zdens(nmats)           ! P'(t), alias for g_zdens_mo 
      integer g_zdens_pdt(nmats)       ! P'(t+dt)
      integer g_zdens_pdt_old(nmats)   ! P'(t+dt) from previous iteration
      integer g_zfock_p12dt(nmats)     ! F'(t+1/2*dt); first extrapolated, then interpolated)
      integer g_zfock_pdt(nmats)       ! F'(t+dt); built from P'(t+dt)
      integer g_zfock_ao(nmats), g_zdens_ao(nmats)
      integer g_fockre(nmats), g_fockim(nmats)
      integer g_movecs_gs(nmats)

      double precision norm1(nmats), norm2(nmats), diff(nmats)
      double complex zjunk

      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid ()

      
C     
C     Determine number of matricies and check consistency.
C
      if ( (nmats .ne. 1).and.(nmats .ne. 2) )
     $     call errquit (pname//"Only works for 1 or 2 matrices",0,0)

      do imat = 1, nmats
       write (amat(imat), "(i1)") imat !assumes single digit number of mats (ok since max 2)
      enddo

      
C
C     Allocation
C
      ok = .true.

      do imat = 1, nmats

         ok = ok .and. ga_create (mt_dcpl, params%ns_ao, params%ns_ao,
     $        "F in magnus inter "//amat(imat), 0, 0, g_zfock_ao(imat))

         ok = ok .and. ga_create (mt_dbl, params%ns_mo, params%ns_mo,
     $        "real fock "//amat(imat), 0, 0, g_fockre(imat))

         ok = ok .and. ga_create (mt_dbl, params%ns_mo, params%ns_mo,
     $        "imaginary fock "//amat(imat), 0, 0, g_fockim(imat))

         ok = ok .and. ga_create (mt_dbl, params%ns_mo, params%ns_mo,
     $        "vecs gs "//amat(imat), 0, 0, g_movecs_gs(imat))

      end do
      if (.not. ok) call errquit (pname//"allocation failed",0,GA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,params%ns_mo,'alphaeval',levals(1),
     &  ievals(1)))
     &   call errquit('cheb_poly: cannot allocate alphaeval',0, MA_ERR)

      if (nmats.eq.2) then
       if (.not.MA_Push_Get(MT_Dbl,params%ns_mo,'betaeval',levals(2),
     &  ievals(2)))
     &     call errquit('cheb_poly: cannot allocate betaeval',0, MA_ERR)
       end if             
 

C
C
C     (2) Propagate P'(t) -> P'(t+dt)
C
         ok = tdfock (params, tt, g_zdens_ao, energies, g_zfock_ao)
         
         if (.not.ok) call errquit (pname//"tdfock failed", 0, 0)
         do imat = 1, nmats
            call canorg_trans (params, "F", "AO->MO",
     $           g_zfock_ao(imat), g_zfock(imat))
         end do

         if (it.eq.1) then
          do imat = 1, nmats
           call convert_z2d (g_zfock(imat), g_fockre(imat),
     &      g_fockim(imat))
           call ga_diag_std (g_fockre(imat), g_movecs_gs(imat),
     &      dbl_mb(ievals(imat)))
           scalefock(imat) = dbl_mb(ievals(imat)+params%ns_ao)
     &     -dbl_mb(ievals(imat)+0)   
c        write(6,*)"scale",scalefock(imat)  
           end do 
          end if
 
        do imat = 1, nmats
         call prop_cheb(g_zdens_mo(imat),g_zfock(imat),
     &    scalefock(imat),params%ns_mo,dt)
        end do



C
C     Clean up.
C
      if (nmats.eq.2) then 
       if (.not.ma_pop_stack(levals(2)))
     &  call errquit('prop_cheb_poly: cannot pop stack betaevals',3, 
     &  MA_ERR)   
       end if

       if (.not.ma_pop_stack(levals(1)))
     &  call errquit('prop_cheb_poly: cannot pop stack alphaevals',3, 
     &  MA_ERR)   
       
      do imat = 1, nmats          
         ok = .true.
         ok = ok .and. ga_destroy (g_movecs_gs(imat))
         ok = ok .and. ga_destroy (g_fockim(imat))
         ok = ok .and. ga_destroy (g_fockre(imat))
         ok = ok .and. ga_destroy (g_zfock_ao(imat))
      
      if (.not. ok) call errquit (pname//"destroys failed",0,GA_ERR)
       end do
      


C
C     If profiling is enabled print convergence data.
C      
      if (params%prof) then 
         write (outstring,"(a,i0,a)")
     $        "Chebyshev propagation converged after ", iinter,
     $        "terms;"
         call prof_end (elapsed, trim(outstring))
      endif

      end subroutine


         
c $Id$
