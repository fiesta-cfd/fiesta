/*
  Copyright 2019-2021 The University of New Mexico

  This file is part of FIESTA.
  
  FIESTA is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option) any
  later version.
  
  FIESTA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
  details.
  
  You should have received a copy of the GNU Lesser General Public License
  along with FIESTA.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef LSDEBUG_H
#define LSDEBUG_H
#define MYDBG printf("%s:%d\n", __FILE__, __LINE__);
#define MYDBGR printf("%d - %s:%d\n", cf.rank, __FILE__, __LINE__);
#define MYDBG0                                                                 \
  if (cf.rank == 0)                                                            \
    printf("%s:%d\n", __FILE__, __LINE__);
#endif
