find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_library(GMP_LIBRARY NAMES gmp libgmp)
set(GMP_LIBRARIES ${GMP_LIBRARY})
