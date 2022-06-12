#include <vector>

using namespace std;

const int Q=9;

const vector<double> w {4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};

// #define x_dim 256
// #define y_dim 256

class D2Q9_LBM{
    int x_dim;
    int y_dim;
    double U;//顶盖的速度
    double rho0;//密度
    double Re;//雷诺数
    double L;//方腔长度
    double v;//运动粘度
    double tao;//松弛时间
    // double f[x_dim][y_dim][Q]; 
    // double u[x_dim][y_dim][2];
    vector<vector<double>> rho;

    vector<vector<vector<double>>> f0; //平衡分布

    vector<vector<vector<double>>> f;

    vector<vector<vector<double>>> fs; //流动后的分布函数，不能直接在f上应用流动函数 
    //u[i][j][0]->x  u[i][j][1]->y
    vector<vector<vector<double>>> u; //宏观速度

    vector<vector<double>> sf; //流函数值
    
    double e[Q][2]={           //代表九个速度方向
        {0,0},
        {1,0},
        {0,1},
        {-1,0},
        {0,-1},
        {1,1},
        {-1,1},
        {-1,-1},
        {1,-1},
    };


    public:

    D2Q9_LBM(int x_dim,int y_dim,double U,double rho0,double Re,double L)\
        :x_dim(x_dim),y_dim(y_dim),rho0(rho0),U(U),Re(Re),L(L){}

    ~D2Q9_LBM(){}

    //初始化分布函数以及初始速度，并且计算计算松弛时间tao
    void init();

    void stream();

    void boundary_processing();

    void calculate_macro_quantities();

    void collision();

};