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
add_library(sparskit ${sparskit_source})
install(TARGETS sparskit EXPORT lasConfig DESTINATION lib)
