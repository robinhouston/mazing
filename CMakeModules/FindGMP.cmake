find_path(GMP_INCLUDE_DIR NAMES gmp.h DOC "The directory containing the GMP header files")
find_library(GMP_LIBRARY NAMES gmp libgmp DOC "The directory containing the GMP library")
set(GMP_LIBRARIES ${GMP_LIBRARY})
