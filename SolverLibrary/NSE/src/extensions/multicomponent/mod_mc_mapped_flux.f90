module mod_mc_mapped_flux
  use mod_precision, only: dp
  use mod_mc_state_layout, only: mc_state_layout
  use mod_mc_euler_config, only: mc_euler_config
  use mod_mc_geometry
  use mod_mc_boundary, only: mc_boundary_state, mc_boundary_is_periodic
  use mod_mc_euler_field, only: evaluate_mc_primitive
  use mod_mc_mapped_transport_provider, only: transport_flux
  implicit none
  private
  public :: mc_mapped_rhs, mc_mapped_rate
contains
  subroutine primitive_vector(q,l,c,v,rho,p,sound)
    real(dp),intent(in)::q(:)
    type(mc_state_layout),intent(in)::l
    type(mc_euler_config),intent(in)::c
    real(dp),intent(out)::v(:),rho,p,sound
    call evaluate_mc_primitive(q,l,c%gamma,rho,v(1:3),v(4),v(5:),p,sound)
  end subroutine
  subroutine mapped_boundary(q,g,l,c,face,area)
    real(dp),intent(in)::q(:),area(3)
    real(dp),intent(out)::g(:)
    type(mc_state_layout),intent(in)::l
    type(mc_euler_config),intent(in)::c
    integer,intent(in)::face
    type(mc_euler_config)::rotated
    real(dp)::t(3,3),n(3),s(size(q))
    integer::f
    t=0
    t(1,1)=1
    t(2,2)=1
    t(3,3)=1
    if((face+1)/2 == 2) then
      n=area/sqrt(sum(area**2))
      t(1,:)=[n(2),-n(1),0.0_dp]
      t(2,:)=n
    end if
    s=q
    s(l%momentum)=matmul(t,q(l%momentum))
    rotated=c
    do f=1,6
      rotated%boundary_reference_velocities(:,f)=matmul(t,c%boundary_reference_velocities(:,f))
    end do
    call mc_boundary_state(s,g,l,rotated,face)
    g(l%momentum)=matmul(transpose(t),g(l%momentum))
  end subroutine
  subroutine face_states(q,l,c,idx,d,ql,qr,left,right,boundary)
    real(dp),intent(in)::q(:,:,:,:)
    type(mc_state_layout),intent(in)::l
    type(mc_euler_config),intent(in)::c
    integer,intent(in)::idx(3),d
    real(dp),intent(out)::ql(:),qr(:)
    integer,intent(out)::left(3),right(3),boundary
    integer::extent(3)
    real(dp)::area(3)
    extent=shape(q(:,:,:,1))
    left=idx
    right=idx
    right(d)=idx(d)+1
    boundary=0
    area=mc_face_area(c,idx,d)
    if(idx(d)==0) then
      right(d)=1
      qr=q(right(1),right(2),right(3),:)
      if(mc_boundary_is_periodic(c,2*d-1)) then
        left(d)=extent(d)
        ql=q(left(1),left(2),left(3),:)
      else
        left=right
        boundary=2*d-1
        call mapped_boundary(qr,ql,l,c,boundary,area)
      end if
    else if(idx(d)==extent(d)) then
      ql=q(left(1),left(2),left(3),:)
      if(mc_boundary_is_periodic(c,2*d)) then
        right(d)=1
        qr=q(right(1),right(2),right(3),:)
      else
        right=left
        boundary=2*d
        call mapped_boundary(ql,qr,l,c,boundary,area)
      end if
    else
      ql=q(left(1),left(2),left(3),:)
      qr=q(right(1),right(2),right(3),:)
    end if
  end subroutine
  subroutine gradients(q,l,c,g)
    real(dp),intent(in)::q(:,:,:,:)
    type(mc_state_layout),intent(in)::l
    type(mc_euler_config),intent(in)::c
    real(dp),intent(out)::g(:,:,:,:,:)
    integer::i,j,k,d,idx(3),a(3),b(3),wall
    real(dp)::ql(l%nvariables),qr(l%nvariables),vp(l%nspecies+4),vm(l%nspecies+4)
    real(dp)::rho,p,sound,h(3),raw(l%nspecies+4,3),m(3,3)
    h=mc_spacing(c)
    !$omp parallel do collapse(3) default(shared) schedule(static) &
    !$omp private(i,j,k,d,idx,a,b,wall,ql,qr,vp,vm,rho,p,sound,raw,m)
    do k=1,c%nz
      do j=1,c%ny
        do i=1,c%nx
          do d=1,3
            idx=[i,j,k]
            call face_states(q,l,c,idx,d,ql,qr,a,b,wall)
            call primitive_vector(qr,l,c,vp,rho,p,sound)
            idx(d)=idx(d)-1
            call face_states(q,l,c,idx,d,ql,qr,a,b,wall)
            call primitive_vector(ql,l,c,vm,rho,p,sound)
            raw(:,d)=(vp-vm)/(2*h(d))
          end do
          m=mc_metric_inverse(c,i,j,k)
          g(i,j,k,:,:)=matmul(raw,m)
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine
  subroutine mc_mapped_rhs(q,rhs,l,c,transport)
    real(dp),intent(in)::q(:,:,:,:)
    real(dp),intent(out)::rhs(:,:,:,:)
    type(mc_state_layout),intent(in)::l
    type(mc_euler_config),intent(in)::c
    logical,intent(in)::transport
    real(dp),allocatable::grad(:,:,:,:,:)
    real(dp)::ql(l%nvariables),qr(l%nvariables),fl(l%nvariables),fr(l%nvariables),flux(l%nvariables)
    real(dp)::vl(l%nspecies+4),vr(l%nspecies+4),gf(l%nspecies+4,3),correction(l%nspecies+4)
    real(dp)::rl,rr,pl,pr,cl,cr,ul,ur,area(3),normal(3),amag,vol,sign,distance,delta(3)
    integer::d,a,b,ia,ib,n,idx(3),left(3),right(3),extent(3),wall,s
    logical::periodic
    extent=shape(q(:,:,:,1))
    rhs=0
    if(transport) then
      allocate(grad(c%nx,c%ny,c%nz,l%nspecies+4,3))
      call gradients(q,l,c,grad)
    else
      allocate(grad(0,0,0,0,0))
    end if
    do d=1,3
      a=mod(d,3)+1
      b=mod(d+1,3)+1
      periodic=mc_boundary_is_periodic(c,2*d)
      !$omp parallel do collapse(2) default(shared) schedule(static) &
      !$omp private(ia,ib,n,idx,left,right,wall,ql,qr,vl,vr,rl,rr,pl,pr,cl,cr,area,normal,amag, &
      !$omp ul,ur,fl,fr,flux,gf,correction,delta,distance,s,vol,sign)
      do ib=1,extent(b)
        do ia=1,extent(a)
          do n=0,extent(d)
            if(periodic .and. n==extent(d)) cycle
            idx=1
            idx(a)=ia
            idx(b)=ib
            idx(d)=n
            call face_states(q,l,c,idx,d,ql,qr,left,right,wall)
            call primitive_vector(ql,l,c,vl,rl,pl,cl)
            call primitive_vector(qr,l,c,vr,rr,pr,cr)
            area=mc_face_area(c,idx,d)
            amag=sqrt(sum(area**2))
            normal=area/amag
            if(transport) then
              gf=0.5_dp*(grad(left(1),left(2),left(3),:,:)+grad(right(1),right(2),right(3),:,:))
              if(wall/=0) then
                distance=mc_cell_volume(c,left(1),left(2),left(3))/amag
                delta=normal*distance
              else
                delta=mc_cell_center(c,right(1),right(2),right(3))- &
                  mc_cell_center(c,left(1),left(2),left(3))
                if(n==0) then
                  ! Periodic seam: use local normal distance, not the box-length jump.
                  distance=0.5_dp*(mc_cell_volume(c,left(1),left(2),left(3))+ &
                    mc_cell_volume(c,right(1),right(2),right(3)))/amag
                  delta=normal*distance
                end if
              end if
              distance=dot_product(delta,normal)
              if(distance<=0) error stop 'non-positive mapped face distance'
              correction=(vr-vl-matmul(gf,delta))/distance
              do s=1,3
                gf(:,s)=gf(:,s)+correction*normal(s)
              end do
              call transport_flux(0.5_dp*(vl+vr),0.5_dp*(rl+rr),gf,l,normal,flux)
              if(wall/=0) then
                if(c%boundary_face_types(wall)=='reflective') then
                  ! Free-slip, adiabatic, impermeable stationary wall.
                  flux(l%first_species:l%last_species)=0
                  flux(l%momentum)=normal*dot_product(flux(l%momentum),normal)
                  flux(l%total_energy)=0
                end if
              end if
              sign=1
            else
              ul=dot_product(vl(1:3),normal)
              ur=dot_product(vr(1:3),normal)
              fl=ql*ul
              fr=qr*ur
              fl(l%momentum)=fl(l%momentum)+pl*normal
              fr(l%momentum)=fr(l%momentum)+pr*normal
              fl(l%total_energy)=(ql(l%total_energy)+pl)*ul
              fr(l%total_energy)=(qr(l%total_energy)+pr)*ur
              flux=0.5_dp*(fl+fr)-0.5_dp*max(abs(ul)+cl,abs(ur)+cr)*(qr-ql)
              sign=-1
            end if
            flux=sign*amag*flux
            if(wall==0 .or. mod(wall,2)==0) then
              vol=mc_cell_volume(c,left(1),left(2),left(3))
              rhs(left(1),left(2),left(3),:)=rhs(left(1),left(2),left(3),:)+flux/vol
            end if
            if(wall==0 .or. mod(wall,2)==1) then
              vol=mc_cell_volume(c,right(1),right(2),right(3))
              rhs(right(1),right(2),right(3),:)=rhs(right(1),right(2),right(3),:)-flux/vol
            end if
          end do
        end do
      end do
      !$omp end parallel do
    end do
  end subroutine
  real(dp) function mc_mapped_rate(c,i,j,k,u,sound) result(rate)
    type(mc_euler_config),intent(in)::c
    integer,intent(in)::i,j,k
    real(dp),intent(in)::u(3),sound
    integer::d,side,p(3)
    real(dp)::a(3)
    rate=0
    do d=1,3
      do side=0,1
        p=[i,j,k]
        p(d)=p(d)-side
        a=mc_face_area(c,p,d)
        rate=rate+0.5_dp*(abs(dot_product(u,a))+sound*sqrt(sum(a*a)))
      end do
    end do
    rate=rate/mc_cell_volume(c,i,j,k)
  end function
end module mod_mc_mapped_flux
