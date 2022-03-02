set(sparskit_source
  src/blas1.f
  src/ilu.f
  src/itaux.f
  src/iters.f
  src/ludcmp.f
  src/matvec.f
  src/pgmres.f
  src/pythag.f
  src/tql2.f
  src/tred2EISPACK.f
  src/utils.f )
add_library(sparskit STATIC ${sparskit_source})
target_compile_options(sparskit PUBLIC
  $<$<AND:$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>,$<CONFIG:DEBUG>>:-Wall -Wextra>
  $<$<AND:$<COMPILE_LANG_AND_ID:C,GNU,Clang>,$<CONFIG:DEBUG>>:-Wall -Wextra>
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU,Clang>:-fdefault-real-8 -fdefault-double-8 -Wall -cpp -std=legacy>
  $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-fdefault-real-8 -fdefault-double-8 -Wall -cpp>)
set_property(TARGET sparskit PROPERTY POSITION_INDEPENDENT_CODE ON)

install(TARGETS sparskit EXPORT las_las DESTINATION lib)
