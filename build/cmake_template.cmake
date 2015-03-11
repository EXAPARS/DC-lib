#  Copyright 2014 - UVSQ
#  Authors list: Loïc Thébault, Eric Petit
#
#  This file is part of the DC-lib.
#
#  DC-lib is free software: you can redistribute it and/or modify it under the
#  terms of the GNU Lesser General Public License as published by the Free Software
#  Foundation, either version 3 of the License, or (at your option) any later version.
#
#  DC-lib is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with
#  the DC-lib. If not, see <http://www.gnu.org/licenses/>.

SET(CMAKE_SYSTEM_NAME Linux)
SET(CMAKE_SYSTEM_PROCESSOR k1om)
SET(CMAKE_SYSTEM_VERSION 1)

# specify the cross compiler
SET(CMAKE_C_COMPILER   icc)
SET(CMAKE_CXX_COMPILER icpc)
SET(MPI_C_COMPILER mpiicc)
SET(_CMAKE_TOOLCHAIN_PREFIX  x86_64-k1om-linux-)

# where is the target environment 
SET(CMAKE_FIND_ROOT_PATH /usr/linux-k1om-4.7)

