struct computeMetrics2D {

    FS4D metrics;
    FS4D grid;

    computeMetrics2D (FS4D m_, FS4D g_) : metrics(m_), grid(g_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        double x_xi, y_xi, x_et, y_et, jac;
        int ii, jj;

        ii = i-3;
        jj = j-3;

       x_xi = (grid(ii+1,jj  ,0,0)+grid(ii+1,jj+1,0,0))/2.0
             -(grid(ii  ,jj  ,0,0)+grid(ii  ,jj+1,0,0))/2.0;
       y_xi = (grid(ii+1,jj  ,0,1)+grid(ii+1,jj+1,0,1))/2.0
             -(grid(ii  ,jj  ,0,1)+grid(ii  ,jj+1,0,1))/2.0;

       x_et = (grid(ii  ,jj+1,0,0)+grid(ii+1,jj+1,0,0))/2.0
             -(grid(ii  ,jj  ,0,0)+grid(ii+1,jj  ,0,0))/2.0;
       y_et = (grid(ii  ,jj+1,0,1)+grid(ii+1,jj+1,0,1))/2.0
             -(grid(ii  ,jj  ,0,1)+grid(ii+1,jj  ,0,1))/2.0;

       jac = 1.0/(x_xi*y_et - x_et*y_xi);

       metrics(i,j,0,0) =  jac*y_et;
       metrics(i,j,0,1) = -jac*x_et;
       metrics(i,j,1,0) = -jac*y_xi;
       metrics(i,j,1,1) =  jac*x_xi;

    }
};
