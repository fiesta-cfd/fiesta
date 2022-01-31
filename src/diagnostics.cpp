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

#include "diagnostics.hpp"
#include "kokkosTypes.hpp"
#include "log2.hpp"
#include <map>

std::string vform(std::vector<FSCAL> vec){
  int len = vec.size();
  std::string str="[";
  for (int i=0; i<len-1; ++i){
    str = fmt::format("{}{:.2e},",str,vec[i]);
  }
  str = fmt::format("{}{:.2e}]",str,vec[len-1]);
  return str;
}

Diagnostics::Diagnostics(){}

Diagnostics::Diagnostics(size_t g_, size_t i_, size_t j_, size_t k_, size_t v_, int freq_):ng(g_),ni(i_),nj(j_),nk(k_),nv(v_),freq(freq_){
  diag = FS4D("diag", ni, nj, nk, nv);
}

Diagnostics::Diagnostics(Diagnostics&& d1):ng(d1.ng),ni(d1.ni),nj(d1.nj),nk(d1.nk),nv(d1.nv),freq(d1.freq){
  diag = FS4D("diag", ni, nj, nk, nv);
}

Diagnostics& Diagnostics::operator=(Diagnostics&& d1){
  if (&d1 == this)
    return *this;
  ng = d1.ng;
  ni = d1.ni;
  nj = d1.nj;
  nk = d1.nk;
  nv = d1.nv;
  freq = d1.freq;
  diag = FS4D("diag", ni, nj, nk, nv);
  return *this;
}

Diagnostics& Diagnostics::operator=(const Diagnostics& d1){
  if (&d1 == this)
    return *this;
  ng = d1.ng;
  ni = d1.ni;
  nj = d1.nj;
  nk = d1.nk;
  nv = d1.nv;
  freq = d1.freq;
  diag = FS4D("diag", ni, nj, nk, nv);
  return *this;
}

void Diagnostics::init(size_t t, FS4D &dvar){
  policy_f3 cell_pol = policy_f3( {0,0,0}, {ni, nj, nk});
  for(size_t v=0; v<nv; ++v){
    Kokkos::parallel_for(cell_pol, KOKKOS_LAMBDA(const int i, const int j, const int k){
      diag(i,j,k,v)=0.0;
      dvar(i,j,k,v)=0.0;
    });
  }
}

void Diagnostics::start(size_t t, FS4D &dvar){
  policy_f3 cell_pol = policy_f3( {0,0,0}, {ni, nj, nk});
  for(size_t v=0; v<nv; ++v){
    Kokkos::parallel_for(cell_pol, KOKKOS_LAMBDA(const int i, const int j, const int k){
      diag(i,j,k,v)=dvar(i,j,k,v);
      dvar(i,j,k,v)=0.0;
    });
  }
}

void Diagnostics::stop(std::string name, size_t t, FS4D &dvar, std::map<std::string, FS4D> &dgvar){
  policy_f3 cell_pol = policy_f3( {0,0,0}, {ni, nj, nk});

  FSCAL localmax = 0.0;
  std::vector<FSCAL> max(nv,0);

  if (dgvar.count(name)==0) {
    Log::debug("Creating new diagnostic entry for : '{}'",name);
    dgvar.emplace(name,FS4D(name, ni, nj, nk, nv));
  }

  if (t % freq == 0){
    for(int v=0; v<nv; ++v){
      Kokkos::parallel_reduce(cell_pol, KOKKOS_LAMBDA(const int i, const int j, const int k, FSCAL &m){
        if (dvar(i,j,k,v) > m) m = dvar(i,j,k,v);
      }, Kokkos::Max<FSCAL>(localmax));
      max[v] = localmax;
    }
    Log::info("'{}' {}",name,vform(max));
  }

  FS4D tmp = dgvar[name];
  for(int v=0; v<nv; ++v){
    Kokkos::parallel_for(cell_pol, KOKKOS_LAMBDA(const int i, const int j, const int k){
      tmp(i,j,k,v) = dvar(i,j,k,v);
      dvar(i,j,k,v) = diag(i,j,k,v) + dvar(i,j,k,v);
    });
  }
}

void Diagnostics::finalize(size_t t, FS4D &dvar, std::map<std::string,FS4D> &dgvar){
  policy_f3 cell_pol = policy_f3( {0,0,0}, {ni, nj, nk});

  std::string name{"dvar"};

  if (dgvar.count(name)==0) {
    Log::debug("Creating new diagnostic entry for : '{}'",name);
    dgvar.emplace(name,FS4D(name, ni, nj, nk, nv));
  }

  FS4D tmp = dgvar[name];

  for(int v=0; v<nv; ++v){
    Kokkos::parallel_for(cell_pol, KOKKOS_LAMBDA(const int i, const int j, const int k){
      tmp(i,j,k,v) = dvar(i,j,k,v);
    });
  }

  if (t % freq == 0){
    FSCAL localmax = 0;
    std::vector<FSCAL> max(nv,0);
    for(size_t v=0; v<nv; ++v){
      Kokkos::parallel_reduce(cell_pol, KOKKOS_LAMBDA(const int i, const int j, const int k, FSCAL &m){
        if (dvar(i,j,k,v) > m) m = dvar(i,j,k,v);
      }, Kokkos::Max<FSCAL>(localmax));
      max[v] = localmax;
    }
    Log::info("'{}' {}",name,vform(max));
  }
}
