      subroutine init_dynamics_grid (grid,imglobal,jmglobal,lm,nt,ak,bk)
      use g3_dynamics_state_module
      type ( dynamics_grid_type ) grid

      integer imglobal,jmglobal,lm,nt
      real(kind=8) ak(lm+1), bk(lm+1)

      integer im,jm,L
      real(kind=8) pi,dl,dp

      im = grid%lattice%im( grid%lattice%pei )
      jm = grid%lattice%jm( grid%lattice%pej )

      pi = 4.0*atan(1.0)
      dl = 2.0*pi/ imglobal
      dp =     pi/(jmglobal-1)

      grid%n  = 1
      grid%im = im
      grid%jm = jm
      grid%lm = lm
      grid%ptracer = nt
      grid%ntracer = nt

      grid%lam_np =  0.0
      grid%phi_np = 90.0
      grid%lam_0  =  0.0

      grid%dl      = dl
      grid%dp      = dp
      grid%dlam(:) = dl
      grid%dphi(:) = dp

      grid%ptop    = ak(1)
      grid%alf (:) = ak(:)
      grid%bet (:) = bk(:)
      grid%sige(:) = grid%bet(:)
      do   L=1,lm
      grid%dsig(L) =  grid%sige(L+1) - grid%sige(L)
      grid% sig(L) = (grid%sige(L+1) + grid%sige(L))*0.5
      enddo

      return
      end

