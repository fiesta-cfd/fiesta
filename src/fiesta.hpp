#include "Kokkos_Core.hpp"
//typedef typename Kokkos::LayoutLeft FS_LAYOUT
//typedef typename Kokkos::View<double******, FS_LAYOUT> FS6D;
//typedef typename Kokkos::View<double*****, FS_LAYOUT> FS5D;
//typedef typename Kokkos::View<double****, FS_LAYOUT> FS4D;
//typedef typename Kokkos::View<double***, FS_LAYOUT> FS3D;
//typedef typename Kokkos::View<double**, FS_LAYOUT> FS2D;
//typedef typename Kokkos::View<double*, FS_LAYOUT> FS1D;

//namespace Fiesta {
    typedef typename Kokkos::View<double******> FS6D;
    typedef typename Kokkos::View<double*****> FS5D;
    typedef typename Kokkos::View<double****> FS4D;
    typedef typename Kokkos::View<double***> FS3D;
    typedef typename Kokkos::View<double**> FS2D;
    typedef typename Kokkos::View<double*> FS1D;
//}
