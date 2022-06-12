#include "lbm.hpp"
#include <math.h>
void D2Q9_LBM::init(){

    f0=vector<vector<vector<double>>>(x_dim,vector<vector<double>>(y_dim,vector<double>(Q,0)));
    f=vector<vector<vector<double>>>(x_dim,vector<vector<double>>(y_dim,vector<double>(Q,0)));
    fs=vector<vector<vector<double>>>(x_dim,vector<vector<double>>(y_dim,vector<double>(Q,0)));
    u=vector<vector<vector<double>>>(x_dim,vector<vector<double>>(y_dim,vector<double>(2,0)));
    rho=vector<vector<double>>(x_dim,vector<double>(y_dim));
    sf=vector<vector<double>>(x_dim,vector<double>(y_dim));

    v=U*L/Re;
    tao=(6*v+1.0)/2.0;

    //初始化速度以及密度
    for (int i=0;i<x_dim;i++){
        for (int j=0;j<y_dim;j++){
            // u[i][j][0]=0;
            // u[i][j][1]=0;
            rho[i][j]=rho0;
        }
        u[i][y_dim-1][0]=U;
    }

    //初始化分布函数
    for (int i=0;i<x_dim;i++){
        for (int j=0;j<y_dim;j++){
            for (int k=0;k<Q;k++){
                double f_item=e[k][0]*u[i][j][0]+e[k][1]*u[i][j][1];//e1i*u
                double s_item=9.0/2*pow(f_item,2);
                double t_item=3.0/2*(pow(u[i][j][0],2)+pow(u[i][j][1],2));
                f0[i][j][k]=w[k]*rho0*(1+3*f_item+s_item-t_item);
                f[i][j][k]=f0[i][j][k];
            }
        }
    }
}

//一个时间步的stream，运行一次函数代表一个时间步
void D2Q9_LBM::stream(){
    for (int i=0;i<x_dim;i++){
        for (int j=0;j<y_dim;j++){
            for (int k=0;k<Q;k++){
                int x_idx=i+e[k][0];
                int y_idx=j+e[k][1];
                if (x_idx>=0 && x_idx<x_dim && y_idx>=0 && y_idx<y_dim)
                    fs[x_idx][y_idx][k]=f[i][j][k];
            }
        }
    }
}

//边界处理
void D2Q9_LBM::boundary_processing(){
    //上下边界
    for (int i=1;i<x_dim-1;i++){
        //下
        f[i][0][2]=fs[i][0][4];
        f[i][0][5]=fs[i][0][7];
        f[i][0][6]=fs[i][0][8];

        // f[i][0][4]=0;
        // f[i][0][7]=0;
        // f[i][0][8]=0;


        //上
        for (int k=0;k<Q;k++){
            f[i][y_dim-1][k]=f0[i][y_dim-1][k];
        }
    }

    for(int j=1;j<y_dim-1;j++){
        //左边界
        f[0][j][1]=fs[0][j][3];
        f[0][j][5]=fs[0][j][7];
        f[0][j][8]=fs[0][j][6];

        // f[0][j][3]=0;
        // f[0][j][6]=0;
        // f[0][j][7]=0;
        
        //右边界
        f[x_dim-1][j][6]=fs[x_dim-1][j][8];
        f[x_dim-1][j][3]=fs[x_dim-1][j][1];
        f[x_dim-1][j][7]=fs[x_dim-1][j][5];

        // f[x_dim-1][j][8]=0;
        // f[x_dim-1][j][1]=0;
        // f[x_dim-1][j][5]=0;
    }

    //四角
    f[0][0][5]=fs[0][0][7];
    f[x_dim-1][y_dim-1][7]=fs[x_dim-1][y_dim-1][5];
    f[0][y_dim-1][8]=fs[0][y_dim-1][6];
    f[x_dim-1][0][6]=fs[x_dim-1][0][8];

}

void D2Q9_LBM::calculate_macro_quantities(){
    for (int i=0;i<x_dim;i++){
        for(int j=0;j<y_dim;j++){
            rho[i][j]=0;
            u[i][j][0]=0;
            u[i][j][1]=0;
            for(int k=0;k<Q;k++){
                rho[i][j]+=f[j][j][k];
                u[i][j][0]+=f[i][j][k]*e[k][0];
                u[i][j][1]+=f[i][j][k]*e[k][1];
            }
        }
    }
}

void D2Q9_LBM::collision(){
    for (int i=0;i<x_dim;i++){
        for (int j=0;j<y_dim;j++){
            for(int k=0;k<Q;k++){
                f[i][j][k]=fs[i][j][k]-1/tao*(f[i][j][k]-f0[i][j][k]);
            }
        }
    }
}