#!/bin/sh

#=======================================================================
#   Copyright (C) 2024 Univ. of Bham  All rights reserved.
#   
#   		FileName：		autogen.sh
#   	 	Author：		LongLI <long.l@cern.ch>
#   		Time：			2024.04.21
#   		Description：
#
#======================================================================

srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.

(cd $srcdir; aclocal -I ${OFFLINE_MAIN}/share;\
 libtoolize --force; automake -a --add-missing; autoconf)

$srcdir/configure "$@"
