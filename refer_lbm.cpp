#include<iostream> 
#include<cmath> 
#include<cstdlib> 
#include<iomanip> 
#include<fstream> 
#include<sstream> 
#include<string> 

#define Q 9       //9个速度方向
#define NX 256   //格子横向个数
#define NY 256   //格子纵向个数
#define rho0 2.7		//初始密度
#define U 0.1			//顶部速度
#define Re 7500      //Re数

using namespace std;

int e[Q][2] = { {0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1} };    //方格子的速度（包括方向）
double w[Q] = { 4.0 / 9,1.0 / 9,1.0 / 9,1.0 / 9,1.0 / 9,1.0 / 36,1.0 / 36,1.0 / 36,1.0 / 36 };      //权系数
double rho[NX + 1][NY + 1];  //位置（X,Y）的密度
double u[NX + 1][NY + 1][2];	//位置（X,Y）的速度（二维）
double u0[NX + 1][NY + 1][2];		//中间变量，用于存储演化后的速度
double f[NX + 1][NY + 1][Q];		//分布函数 f[x坐标][y坐标][速度方向]
double fs[NX + 1][NY + 1][Q];       //流动后的分布函数
double sf[NX + 1][NY + 1]; //流函数值
double maxsf; //10000时间步之后的绝对值最大的流函数值
double maxsf1; //10000时间步之前的绝对值最大的流函数值
double maxsfchazhi;  //连续10000步的流函数差值，用来判断收敛
int i, j, k, ip, jp, n;  //各种参数
double c, dx, dy, Lx, Ly, dt, tau_f, niu, error;    //全局变量

double feq(int k, double rho, double u[2]);    //均衡函数
void init();        //初始化函数
void stream();          //流动函数
void calc();        //计算宏观量函数
void collide();         //碰撞函数
void boundary();        //边界处理函数
void output2(int m);    //输出速度型数据
void output(int m);     //输出tecplot数据
void Error();       //误差函数
void str_funct(); //流函数


int main()
{
    dx = 1.0;      //x方向的网格步长
    dy = 1.0;			//y方向的网格步长
    Lx = dx * double(NY);   //x方向的长度
    Ly = dy * double(NX);    //y方向的长度
    dt = dx;     //时间步长
    c = dx / dt;     //格子速度	
    niu = U * Lx / Re;    //U是顶盖速度，niu是运动粘度系数
    tau_f = 3.0 * niu + 0.5;   //无量纲松弛时间


    init();//初始化

    for (n = 0; ; n++)
    {
        stream();  //迁移
        boundary(); //边界处理
        calc();  //计算宏观量
        collide();  //碰撞
      //  if(n==0)  //计算初始流场流函数最大值
        str_funct();
        if (n == 0)
        {
            str_funct();
            maxsf1 = maxsf;
        }

        if (n % 100 == 0)   //每100时间步刷新显示计算结果
        {
            Error();
            cout << "第" << n << "次计算的结果:" << endl << "NX/2,NY/2点处的u,v 速度分别为:" << setprecision(7) << u[NX / 2][NY / 2][0] << "," << u[NX / 2][NY / 2][1] << endl;
            cout << "uv的误差为:" << setiosflags(ios::scientific) << error << endl;
            cout << "第" << n << "次计算的流函数绝对值最大的流函数值为：" << maxsf << endl;
            if (n >= 1000)
            {
                if (n % 10000 == 0)    //每10000步输出一次tecplot数据
                {
                    output(n);   //输出tecplot文件
                }



                if (n % 10000 == 0)   //每10000步计算一次流函数绝对值的最大值，与前10000次做差，用来判断收敛
                {

                    str_funct();
                    maxsfchazhi = abs(maxsf1 - maxsf);
                    maxsf1 = maxsf;
                    if (maxsfchazhi < 10e-5)  //如果满足收敛条件，输出速度型文件，结束计算
                    {
                        cout << "计算完成！！！" << endl;
                        cout << "此时连续10000步的流函数差值为" << maxsfchazhi << endl;
                        output2(n);
                        break;
                    }
                }
                if (n == 1100000)    //如果计算1100000次都不收敛，视为不稳定，结束计算
                {
                    cout << "计算可能达不到稳态，已经终止计算！" << endl;
                    break;
                }



            }
        }


    }

}

double feq(int k, double rho, double u[2])     //均衡函数
{
    double eu, uv, feq;
    eu = e[k][0] * u[0] + e[k][1] * u[1];
    uv = u[0] * u[0] + u[1] * u[1];
    feq = w[k] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
    return feq;
}

void init()    //初始化函数
{
    for (i = 0; i <= NX; i++)  //初始化 
        for (j = 0; j <= NY; j++)
        {
            u[i][j][0] = 0;
            u[i][j][1] = 0;
            rho[i][j] = rho0;
            u[i][NY][0] = U;
            for (k = 0; k < Q; k++)
            {
                f[i][j][k] = feq(k, rho[i][j], u[i][j]);
            }
        }
}

void stream()  //迁移（流动）函数    节点中心向四周迁移
{
    for (i = 0; i <= NX; i++)
        for (j = 0; j <= NY; j++)
            for (k = 0; k < Q; k++)
            {
                ip = i + e[k][0];
                jp = j + e[k][1];
                if (ip >= 0 && ip <= NX && jp >= 0 && jp <= NX)    //让边界的可以迁移的方向也参与迁移
                    fs[ip][jp][k] = f[i][j][k];   //fs[ip][jp][k]为迁移后的分布函数
            }

}

void calc()   //计算宏观量
{
    for (i = 1; i < NX; i++)    //计算宏观量,速度密度等 
        for (j = 1; j < NY; j++)
        {
            u0[i][j][0] = u[i][j][0];
            u0[i][j][1] = u[i][j][1];
            rho[i][j] = 0;
            u[i][j][0] = 0;
            u[i][j][1] = 0;
            for (k = 0; k < Q; k++)
            {
                rho[i][j] += fs[i][j][k];
                u[i][j][0] += e[k][0] * fs[i][j][k];
                u[i][j][1] += e[k][1] * fs[i][j][k];
            }
            u[i][j][0] /= rho[i][j];
            u[i][j][1] /= rho[i][j];
        }

}

void collide()  //碰撞函数
{
    for (i = 1; i < NX; i++)
        for (j = 1; j < NY; j++)
            for (k = 0; k < Q; k++)
            {
                f[i][j][k] = fs[i][j][k] + (feq(k, rho[i][j], u[i][j]) - fs[i][j][k]) / tau_f;
            }
}

void boundary()    //边界处理，顶盖均衡，碰墙反弹
{
    for (j = 1; j < NY; j++)
    {

        f[NX][j][3] = fs[NX][j][1];  //右边界          
        f[NX][j][6] = fs[NX][j][8];
        f[NX][j][7] = fs[NX][j][5];

        f[0][j][1] = fs[0][j][3];    //左边界
        f[0][j][5] = fs[0][j][7];
        f[0][j][8] = fs[0][j][6];

    }

    for (i = 1; i < NX; i++)  //上下边界 
        for (k = 0; k < Q; k++)
        {
            u[i][NY][0] = U;    //上边界分布函数设置为平衡分布
            u[i][NY][1] = 0;
            rho[i][NY] = rho0;
            f[i][NY][k] = feq(k, rho[i][NY], u[i][NY]);


            f[i][0][2] = fs[i][0][4];  //下边界
            f[i][0][5] = fs[i][0][7];
            f[i][0][6] = fs[i][0][8];

        }
    f[0][0][5] = fs[0][0][7];       //四个角流入流场的分布函数
    f[NX][0][6] = fs[NX][0][8];
    f[0][NY][8] = fs[0][NY][6];
    f[NX][NY][7] = fs[NX][NY][5];

}

void output(int m) //输出tecplot文件
{
    ostringstream name;
    name << "the" << m << "th" << Re << ").dat";
    ofstream out(name.str().c_str());
    out << "Title=\"LBM Lid Driven Flow\"\n" << "VARIABLES = \"X\",\"Y\",\"U\",\"V\",\"S\"\n" << "ZONE T=\"BOX\",I=" << NX + 1 << ",J=" << NY + 1 << ",F=POINT" << endl;
    for (j = 0; j <= NY; j++)
        for (i = 0; i <= NX; i++)
        {
            out << double(i) / NX << " " << double(j) / NY << " " << u[i][j][0] << " " << u[i][j][1] << " " << sf[i][j] << endl;
        }
}

void output2(int m) //输出速度型TXT文件
{
    ostringstream name;
    name << Re << "vv.txt";
    ofstream out(name.str().c_str());
    for (j = 0; j <= NY; j++)
        out << double(j) / NY << " " << u[j][NX / 2][1] / U << endl;   //输出V的速度型

    ostringstream name2;
    name2 << Re << "uu.txt";
    ofstream out2(name2.str().c_str());
    for (j = 0; j <= NY; j++)
        out2 << double(j) / NY << " " << u[NX / 2][j][0] / U << endl;   //输出U的速度型
}


void Error()   //速度误差函数
{
    double temp1, temp2;
    temp1 = 0;
    temp2 = 0;
    for (i = 1; i < NX; i++)
        for (j = 1; j < NY; j++)
        {
            temp1 += ((u[i][j][0] - u0[i][j][0]) * (u[i][j][0] - u0[i][j][0]) + (u[i][j][1] - u0[i][j][1]) * (u[i][j][1] - u0[i][j][1]));
            temp2 += (u[i][j][0] * u[i][j][0] + u[i][j][1] * u[i][j][1]);
        }
    temp1 = sqrt(temp1);
    temp2 = sqrt(temp2);
    error = temp1 / temp2;
}

void str_funct()   //计算流场流函数的绝对值最大的点
{
    for (j = 0; j <= NY; j++)
        for (i = 1; i <= NX; i++)
        {
            if (i == 1) sf[1][j] = -(u[0][j][1] + u[1][j][1]) / 2.0 * dx / U / Lx;
            if (i >= 2) sf[i][j] = sf[i - 2][j] - (u[i - 2][j][1] + 4 * u[i - 1][j][1] + u[i][j][1]) / 3.0 * dx / U / Lx;
        }
    maxsf = sf[1][1];
    for (j = 1; j < NY; j++)
        for (i = 1; i < NX; i++)
        {
            if (abs(sf[j][i]) > abs(maxsf))
                maxsf = sf[j][i];
        }
}

