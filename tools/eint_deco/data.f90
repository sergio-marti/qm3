program const
    use dynamo
    implicit none
    call dynamo_header
    call mm_file_process( "borra", "opls_amber" )
    call mm_system_construct( "borra", "seq" )
    open( unit = 666, file = "data", action = "write", form = "unformatted" )
    ! mm_system.F90 -------------------------------------------
    ! ATMEPS(I)   = 2.0_DP * SQRT ( TYPES(MM)%EPSILON )
    write( 666 ) atmeps * 0.5d0
    ! ATMSIG(I)   = 0.5_DP * TYPES(MM)%SIGMA
    write( 666 ) ( 2.0d0 * atmsig ) * 0.5612310241546865d0
    ! ---------------------------------------------------------
    close( 666 )
end
