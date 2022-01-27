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

#include "kokkosTypes.hpp"
#include "input.hpp"
#include "diagnostics.hpp"
#include <map>

void diagnostics::start(fsconf cf,FS4D dvar){
  policy_f3 cell_pol = policy_f3( {cf.ng, cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk - cf.ng});
  for(int v=0; v<cf.nvt; ++v){
    Kokkos::parallel_for(cell_pol, KOKKOS_LAMBDA(const int i, const int j, const int k){
      dvar(i,j,k,v)=0.0;
      diag(i,j,k,v)=0.0;
    });
  }
}

void diagnostics::check(std::string name, FS4D dvar, std::map<std::string, FS4D> dgvar){
  if (cf.t % freq == 0){
    policy_f3 cell_pol = policy_f3( {cf.ng, cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk - cf.ng});

    FSCAL localmax;
    std::vector<FSCAL> max(cf.nvt,0);

    for(int v=0; v<cf.nvt; ++v){
      Kokkos::parallel_reduce(cell_pol, KOKKOS_LAMBDA(const int i, const int j, const int k, FSCAL &m){
        if (dvar(i,j,k,v) > m) m = dvar(i,j,k,v);
      }, Kokkos::Max<FSCAL>(localmax));
      #ifdef HAVE_MPI
        MPI_Allreduce(&localmax, &max[v], 1, MPI_DOUBLE, MPI_MAX, cf.comm);
      #else
        max[v] = localmax;
      #endif
    }
    Log::infoWarning("'{}' {}",name,max);

    for(int v=0; v<cf.nvt; ++v){
      Kokkos::parallel_for(cell_pol, KOKKOS_LAMBDA(const int i, const int j, const int k){
        diag(i,j,k,v) += dvar(i,j,k,v);
        dvar(i,j,k,v) = 0.0;
      });
    }
  }
}

void diagnostics::stop(fsconf cf,FS4D dvar){
  policy_f3 cell_pol = policy_f3( {cf.ng, cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk - cf.ng});
  for(int v=0; v<cf.nvt; ++v){
    Kokkos::parallel_for(cell_pol, KOKKOS_LAMBDA(const int i, const int j, const int k){
      dvar(i,j,k,v) = diag(i,j,k,v);
    });
  }
}
