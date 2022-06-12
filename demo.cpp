#include"lbm.hpp"

using namespace std;

int main(){
    double U=0.1;
    double rho0=2.7;
    int x_dim=256;
    int y_dim=256;
    double Re=7400;
    double L=256;
    // double v;
    // double tao;

    D2Q9_LBM lbm_demo(x_dim,y_dim,U,rho0,Re,L);
    int n=0;
    lbm_demo.init();
    while(n<10000){
        lbm_demo.stream();
        lbm_demo.boundary_processing();
        lbm_demo.calculate_macro_quantities();
        lbm_demo.collision();
        n++;
    }
}