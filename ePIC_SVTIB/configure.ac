AC_INIT(ePIC_SVTIB, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

case $CXX in
	clang++)
	CXXFLAGS="$CXXFLAGS -Wall -Werror -Wno-undefined-var-template"
	;;
	g++)
	CXXFLAGS="$CXXFLAGS -Wall -Werror"
	;;
esac

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
