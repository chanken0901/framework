program rf_wave_unit
  use mod_rf_wave
  use mod_rf_thermo, only: dp,require
  implicit none
  type(wave_history) :: history
  real(dp) :: p(8),x,v,jump
  logical :: found,valid
  p=1;p(5:)=3
  call sample_wave(p,1._dp,0._dp,0._dp,8._dp,1._dp,history,x,v,jump,found,valid)
  call require(found.and..not.valid.and.x==4.and.jump==2,'Initial wave detection')
  p=1;p(6:)=3
  call sample_wave(p,1._dp,1._dp,0._dp,8._dp,1._dp,history,x,v,jump,found,valid)
  call require(found.and.valid.and.x==5.and.v==1,'Right-moving wave speed')
  p=1
  call sample_wave(p,1._dp,2._dp,0._dp,8._dp,1._dp,history,x,v,jump,found,valid)
  call require(.not.found.and..not.valid,'Uniform state is not a shock')
  p(4:)=3
  call sample_wave(p,1._dp,3._dp,0._dp,8._dp,1._dp,history,x,v,jump,found,valid)
  call require(found.and..not.valid.and.x==3,'Reacquisition resets speed')
  p=1;p(3:)=3
  call sample_wave(p,1._dp,4._dp,0._dp,8._dp,1._dp,history,x,v,jump,found,valid)
  call require(found.and.valid.and.v==-1,'Left-moving wave speed')
  call sample_wave(p,1._dp,5._dp,4._dp,8._dp,1._dp,history,x,v,jump,found,valid)
  call require(.not.found,'Search interval excludes unrelated wave')
  call sample_wave(p,1._dp,6._dp,0._dp,8._dp,3._dp,history,x,v,jump,found,valid)
  call require(.not.found,'Pressure threshold rejects weak feature')
  write(*,*) '[OK] wave position, signed speed, loss and reacquisition'
end program
