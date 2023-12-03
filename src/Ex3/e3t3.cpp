#include "e3t3.h"
#include "ui_e3t3.h"
#include <armadillo>
#include <vector>
using namespace arma;
using namespace std;
#include <QFile>
#include <string>
#ifdef DEBUG
#include <QDebug>
namespace Ex3 {
extern bool Col_Gauss(mat,vec,vector<mat>&,vector<vec>&,const double&);
}
#endif
#ifdef PRINT_LATEX
#include <QDebug>
#endif
#ifdef ZSC_TIME
#include <QDebug>
#include <time.h>
#endif
#include <QMessageBox>
#include <QStandardPaths>
#include <QFileDialog>
#include <QInputDialog>
extern double calStr(string);
namespace _MYFUNCTION {
#ifdef PRINT_LATEX
void print_latex_vec(const vec& x,const unsigned& i)
{
    if(x.empty())
    {
        qDebug()<<i<<"\\\\";
        return;
    }
    string s(to_string(i));
    for(const auto& j:x)(s+='&')+=to_string(j);
    qDebug()<<(s+="\\\\").c_str();
}
#endif
#ifdef DEBUG
void debug_print_vec(const vec& x)
{
    string s;
    for(unsigned i(0);i<x.n_elem;++i)
        (s+=to_string(x.at(i)))+='\t';
    qDebug()<<s.c_str();
}
#endif
void write_vec(QFile& file,const vec& x)
{
    string s;
    for(auto& i:x)s+=(to_string(i))+=',';
    *(s.end()-1)='\n';
    file.write(s.c_str());
}
inline cx_double Sqrt(const double& x)noexcept
{
    return x<0.0?cx_double(0.0,sqrt(-x)):cx_double(sqrt(x),0.0);
}
}
namespace Ex3 {
/*
 * Jacobi迭代法(矩阵法)
 * A:系数矩阵
 * b:常数向量
 * x:迭代初始向量
 * u:迭代过程
 * e:迭代精度
 */
void Jacobi_mat(const mat& A,const vec& b,vec& x,vector<vec>& u,const double& e=1e-6)
{
    vec x1(x),b1(b.n_elem,1,fill::none);
    u.push_back(x);
    mat A1(A.n_rows,A.n_cols,fill::none);
    for(unsigned i(0);i!=x.n_rows;++i)
    {
        if(abs(A.at(i,i))<e)throw "对角线元素为零！";
        A1.at(i,i)=0;
        x.at(i)=b1.at(i)=b.at(i)/A.at(i,i);
        for(unsigned j(0);j!=A.n_cols;++j)
            if(i!=j)
                x.at(i)-=(A1.at(i,j)=A.at(i,j)/A.at(i,i))*x1.at(j);
    }
    u.push_back(x);
    while(norm(x1-x)>e)
    {
        x1=x;
        u.push_back(x=b1-A1*x);
    }
}
/*
 * Jacobi迭代法(分量法)
 * A:系数矩阵
 * b:常数向量
 * x:迭代初始向量
 * u:迭代过程
 * e:迭代精度
 */
void Jacobi_ele(const mat& A,const vec& b,vec& x,vector<vec>& u,const double& e=1e-6)
{
    vec x1(x),b1(b.n_elem,1,fill::none);
    u.push_back(x);
    mat A1(A.n_rows,A.n_cols,fill::none);
    for(unsigned i(0);i!=x.n_rows;++i)
    {
        if(abs(A.at(i,i))<e)throw "对角线元素为零！";
        x.at(i)=b1.at(i)=b.at(i)/A.at(i,i);
        for(unsigned j(0);j!=A.n_cols;++j)
            if(i!=j)
                x.at(i)-=(A1.at(i,j)=A.at(i,j)/A.at(i,i))*x1.at(j);
    }
    u.push_back(x);
    while(norm(x1-x)>e)
    {
        x1=x;
        for(unsigned i(0);i!=x.n_rows;++i)
        {
            x.at(i)=b1.at(i);
            for(unsigned j(0);j!=A.n_cols;++j)
                if(i!=j)
                    x.at(i)-=A1.at(i,j)*x1.at(j);
        }
        u.push_back(x);
    }
}
/*
 * Gauss-Seidel迭代法(矩阵法)
 * A:系数矩阵
 * b:常数向量
 * x:迭代初始向量
 * u:迭代过程
 * e:迭代精度
 */
void Gauss_Seidel_mat(const mat& A,const vec& b,vec& x,vector<vec>& u,const double& e=1e-6)
{
    mat T(inv(trimatl(A))),A1(T*trimatu(A,1));
    vec b1(T*b),x1;
    u.clear();
    u.push_back(x);
    do
    {
        x1=x;
        u.push_back(x=b1-A1*x);
    }while(norm(x-x1)>e);
}
/*
 * Gauss-Seidel迭代法(分量法)
 * A:系数矩阵
 * b:常数向量
 * x:迭代初始向量
 * u:迭代过程
 * e:迭代精度
 */
void Gauss_Seidel_ele(const mat& A,const vec& b,vec& x,vector<vec>& u,const double& e=1e-6)
{
    vec x1(x),b1(b.n_elem,1,fill::none);
    u.push_back(x);
    mat A1(A.n_rows,A.n_cols,fill::none);
    for(unsigned i(0);i!=x.n_rows;++i)
    {
        if(abs(A.at(i,i))<e)throw "对角线元素为零！";
        x.at(i)=b1.at(i)=b.at(i)/A.at(i,i);
        unsigned j(-1);
        while(++j!=i)x.at(i)-=(A1.at(i,j)=A.at(i,j)/A.at(i,i))*x.at(j);
        while(++j!=A.n_cols)x.at(i)-=(A1.at(i,j)=A.at(i,j)/A.at(i,i))*x.at(j);
    }
    u.push_back(x);
    while(norm(x1-x)>e)
    {
        x1=x;
        for(unsigned i(0);i!=x.n_rows;++i)
        {
            x.at(i)=b1.at(i);
            unsigned j(-1);
            while(++j!=i)x.at(i)-=A1.at(i,j)*x.at(j);
            while(++j!=A.n_cols)x.at(i)-=A1.at(i,j)*x.at(j);
        }
        u.push_back(x);
    }
}
/*
 * SOR法(矩阵法)
 * A:系数矩阵
 * b:常数向量
 * w:松弛因子
 * x:迭代初始向量
 * u:迭代过程
 * e:迭代精度
 */
void SOR_mat(const mat& A,const vec& b,const double& w,vec& x,vector<vec>& u,const double& e=1e-6)
{
    mat L(trimatl(A,-1)),D(diagmat(A)),U(trimatu(A,1)),T(inv(trimatl(w*L+D))),A1(T*((1-w)*D-w*U));
    vec b1(T*b*w),x1;
    u.clear();
    u.push_back(x);
    do
    {
        x1=x;
        u.push_back(x=b1+A1*x);
    }while(norm(x-x1)>e);
}
/*
 * SOR法(分量法)
 * A:系数矩阵
 * b:常数向量
 * w:松弛因子
 * x:迭代初始向量
 * u:迭代过程
 * e:迭代精度
 */
void SOR_ele(const mat& A,const vec& b,const double& w,vec& x,vector<vec>& u,const double& e=1e-6)
{
    double w1(1-w);
    vec x1(x),b1(b.n_elem,1,fill::none);
    u.push_back(x);
    mat A1(A.n_rows,A.n_cols,fill::none);
    for(unsigned i(0);i!=x.n_rows;++i)
    {
        if(abs(A.at(i,i))<e)throw "对角线元素为零！";
        (x.at(i)*=w1)+=b1.at(i)=w*b.at(i)/A.at(i,i);
        unsigned j(-1);
        while(++j!=i)x.at(i)-=(A1.at(i,j)=w*A.at(i,j)/A.at(i,i))*x.at(j);
        while(++j!=A.n_cols)x.at(i)-=(A1.at(i,j)=w*A.at(i,j)/A.at(i,i))*x.at(j);
    }
    u.push_back(x);
    while(norm(x1-x)>e)
    {
        x1=x;
        for(unsigned i(0);i!=x.n_rows;++i)
        {
            (x.at(i)*=w1)+=b1.at(i);
            unsigned j(-1);
            while(++j!=i)x.at(i)-=A1.at(i,j)*x.at(j);
            while(++j!=A.n_cols)x.at(i)-=A1.at(i,j)*x.at(j);
        }
        u.push_back(x);
    }
}
/*
 * Jacobi迭代法(分量法)
 * A:系数矩阵
 * b:常数向量
 * x:迭代初始向量
 * p:迭代过程保存路径
 * e:迭代精度
 */
void Jacobi(const mat& A,const vec& b,vec& x,const QString& p,const double& e=1e-6)
{
    QFile file(p);
    if(!file.open(QIODevice::WriteOnly))throw "文件保存失败！";
    vec x1(x),b1(b.n_elem,1,fill::none);
    _MYFUNCTION::write_vec(file,x);
#ifdef PRINT_LATEX
    unsigned k(0);
    _MYFUNCTION::print_latex_vec(x,++k);
#endif
    mat A1(A.n_rows,A.n_cols,fill::none);
    for(unsigned i(0);i!=x.n_rows;++i)
    {
        if(abs(A.at(i,i))<e)throw "对角线元素为零！";
        x.at(i)=b1.at(i)=b.at(i)/A.at(i,i);
        for(unsigned j(0);j!=A.n_cols;++j)
            if(i!=j)
                x.at(i)-=(A1.at(i,j)=A.at(i,j)/A.at(i,i))*x1.at(j);
    }
    _MYFUNCTION::write_vec(file,x);
#ifdef PRINT_LATEX
    _MYFUNCTION::print_latex_vec(x,++k);
#endif
    while(norm(x1-x)>e)
    {
        x1=x;
        for(unsigned i(0);i!=x.n_rows;++i)
        {
            x.at(i)=b1.at(i);
            for(unsigned j(0);j!=A.n_cols;++j)
                if(i!=j)
                    x.at(i)-=A1.at(i,j)*x1.at(j);
        }
        _MYFUNCTION::write_vec(file,x);
#ifdef PRINT_LATEX
        _MYFUNCTION::print_latex_vec(x,++k);
#endif
    }
    file.close();
}
/*
 * Gauss-Seidel迭代法(分量法)
 * A:系数矩阵
 * b:常数向量
 * x:迭代初始向量
 * p:迭代过程保存路径
 * e:迭代精度
 */
void Gauss_Seidel(const mat& A,const vec& b,vec& x,const QString& p,const double& e=1e-6)
{
    QFile file(p);
    if(!file.open(QIODevice::WriteOnly))throw "文件保存失败！";
    vec x1(x),b1(b.n_elem,1,fill::none);
    _MYFUNCTION::write_vec(file,x);
#ifdef PRINT_LATEX
    unsigned k(0);
    _MYFUNCTION::print_latex_vec(x,++k);
#endif
    mat A1(A.n_rows,A.n_cols,fill::none);
    for(unsigned i(0);i!=x.n_rows;++i)
    {
        if(abs(A.at(i,i))<e)throw "对角线元素为零！";
        x.at(i)=b1.at(i)=b.at(i)/A.at(i,i);
        unsigned j(-1);
        while(++j!=i)x.at(i)-=(A1.at(i,j)=A.at(i,j)/A.at(i,i))*x.at(j);
        while(++j!=A.n_cols)x.at(i)-=(A1.at(i,j)=A.at(i,j)/A.at(i,i))*x.at(j);
    }
    _MYFUNCTION::write_vec(file,x);
#ifdef PRINT_LATEX
    _MYFUNCTION::print_latex_vec(x,++k);
#endif
    while(norm(x1-x)>e)
    {
        x1=x;
        for(unsigned i(0);i!=x.n_rows;++i)
        {
            x.at(i)=b1.at(i);
            unsigned j(-1);
            while(++j!=i)x.at(i)-=A1.at(i,j)*x.at(j);
            while(++j!=A.n_cols)x.at(i)-=A1.at(i,j)*x.at(j);
        }
        _MYFUNCTION::write_vec(file,x);
#ifdef PRINT_LATEX
        _MYFUNCTION::print_latex_vec(x,++k);
#endif
    }
    file.close();
}
/*
 * SOR法(分量法)
 * A:系数矩阵
 * b:常数向量
 * w:松弛因子
 * x:迭代初始向量
 * p:迭代过程保存路径
 * e:迭代精度
 */
void SOR(const mat& A,const vec& b,const double& w,vec& x,const QString& p,const double& e=1e-6)
{
    QFile file(p);
    if(!file.open(QIODevice::WriteOnly))throw "文件保存失败！";
    double w1(1-w);
    vec x1(x),b1(b.n_elem,1,fill::none);
    _MYFUNCTION::write_vec(file,x);
#ifdef PRINT_LATEX
    unsigned k(0);
    _MYFUNCTION::print_latex_vec(x,++k);
#endif
    mat A1(A.n_rows,A.n_cols,fill::none);
    for(unsigned i(0);i!=x.n_rows;++i)
    {
        if(abs(A.at(i,i))<e)throw "对角线元素为零！";
        (x.at(i)*=w1)+=b1.at(i)=w*b.at(i)/A.at(i,i);
        unsigned j(-1);
        while(++j!=i)x.at(i)-=(A1.at(i,j)=w*A.at(i,j)/A.at(i,i))*x.at(j);
        while(++j!=A.n_cols)x.at(i)-=(A1.at(i,j)=w*A.at(i,j)/A.at(i,i))*x.at(j);
    }
    _MYFUNCTION::write_vec(file,x);
#ifdef PRINT_LATEX
    _MYFUNCTION::print_latex_vec(x,++k);
#endif
    while(norm(x1-x)>e)
    {
        x1=x;
        for(unsigned i(0);i!=x.n_rows;++i)
        {
            (x.at(i)*=w1)+=b1.at(i);
            unsigned j(-1);
            while(++j!=i)x.at(i)-=A1.at(i,j)*x.at(j);
            while(++j!=A.n_cols)x.at(i)-=A1.at(i,j)*x.at(j);
        }
        _MYFUNCTION::write_vec(file,x);
#ifdef PRINT_LATEX
        _MYFUNCTION::print_latex_vec(x,++k);
#endif
    }
    file.close();
}
/*
 * 最速下降法
 * A:系数矩阵
 * b:常数向量
 * x:迭代初始向量
 * u:迭代过程
 * e:迭代精度
 */
void steepest_descent_method(const mat& A,const mat& b,vec& x,vector<vec>& u,const double& e=1e-8)
{
    u.clear();
    u.push_back(x);
    vec x1;
    do
    {
        vec r(b-A*(x1=x));
        u.push_back(x+=dot(r,r)/dot(r,A*r)*r);
    }while(norm(x1-x)>e);
}
/*
 * 共轭梯度法
 * A:系数矩阵
 * b:常数向量
 * x:迭代初始向量
 * u:迭代过程
 * e:迭代精度
 */
void conjugate_gradient_method(const mat& A,const mat& b,vec& x,vector<vec>& u,const double& e=1e-8)
{
    u.clear();
    u.push_back(x);
    vec x1,r(b-A*x),p(r);
    do
    {
        x1=x;
        mat Ap(A*p);
        double alpha(dot(p,r)/dot(p,Ap));
        u.push_back(x+=alpha*p);
        (p*=-dot(p,A*r)/dot(p,Ap))+=r-=alpha*Ap;
    }while(norm(x1-x)>e);
}
/*
 * 奇异值分解法
 * A:系数矩阵
 * b:常数向量
 * x:初始向量
 * u:迭代过程
 * e:迭代精度
 */
void SVD_method(mat& A,const mat& b,vec& x,vector<vec>& u,const double& e=1e-8)
{
    u.clear();
    mat U,V;
    vec s;
    if(!svd(U,s,V,A,"std"))throw "奇异值分解失败！";
    mat Vt(V.t());
    u.push_back(x);
    for(unsigned i(0);i!=A.n_cols;++i)
        if(abs(s.at(i))<e)
            A-=s.at(i)*U.col(i)*Vt.row(i);
        else
            u.push_back(x+=dot(U.col(i),b)/s.at(i)*V.col(i));
}
/*
 * 预处理共轭梯度法
 * A:系数矩阵
 * b:常数向量
 * x:迭代初始向量
 * u:迭代过程
 * M:迭代方式
 *  0:对角元素法
 *  1:Cholesky分解法
 *  2:对称超松弛迭代法
 * e:迭代精度
 * w:松弛因子
 */
//void pre_conjugate_gradient_method(const mat& A,const mat& b,vec& x,vector<vec>& u,unsigned char M='\0',const double& e=1e-8,const double w=1)
//{
//    cx_mat S;
//    switch (M) {
//    case '\0':
//    {
//        S=diagmat(A);
//        S.for_each(_MYFUNCTION::Sqrt);
//        break;
//    }
//    case '\2':
//    {
//        cx_mat D(diagmat(A/w),zeros(size(A))),U(trimatu(A,1),zeros(size(A)));
//        if(w>2.||w<0.)throw "松弛因子输入有误！";
//        inplace_strans((S=solve(D.transform(_MYFUNCTION::Sqrt),D+U))/=sqrt(2-w));
//        break;
//    }
//    default:
//        if(abs(S.set_size(size(A)).at(0,0)=_MYFUNCTION::Sqrt(A.at(0,0)))<e)throw "除数为零！";
//        for(unsigned i(1);i<A.n_cols;++i)
//            S.at(i,0)=A.at(i,0)/S.at(0,0);
//        for(unsigned k(1);k<A.n_rows;++k)
//        {
//            double h((S.at(k,0)*S.at(k,0)).real());
//            for(unsigned i(1);i!=k;++i)
//                h+=(S.at(k,i)*S.at(k,i)).real();
//            S.at(k,k)=A.at(k,k)<h?A.at(k,k):sqrt(A.at(k,k)-h);
//        }
//    }
//}
}
E3t3::E3t3(QWidget *parent):
    QWidget(parent),
    ui(new Ui::E3t3)
{
    ui->setupUi(this);
#ifdef ZSC_TIME
    const mat A{{8,1,2},{8,7,2},{4,9,9}};
    const vec b{10,18,17},x0(zeros(3));
    const unsigned short end(-1);
    {
        auto start_time=clock();
        for(unsigned short i(0);i!=end;++i)
        {
            vec x(x0);
            vector<vec> u;
            Ex3::Jacobi_ele(A,b,x,u);
        }
        auto end_time=clock();
        qDebug()<<"Jacobi迭代法(分量法)耗时:"<<(end_time-start_time)/65536.0;
    }
    {
        auto start_time=clock();
        for(unsigned short i(0);i!=end;++i)
        {
            vec x(x0);
            vector<vec> u;
            Ex3::Jacobi_mat(A,b,x,u);
        }
        auto end_time=clock();
        qDebug()<<"Jacobi迭代法(矩阵法)耗时:"<<(end_time-start_time)/65536.0;
    }
    {
        auto start_time=clock();
        for(unsigned short i(0);i!=end;++i)
        {
            vec x(x0);
            vector<vec> u;
            Ex3::Gauss_Seidel_ele(A,b,x,u);
        }
        auto end_time=clock();
        qDebug()<<"Gauss_Seidel迭代法(分量法)耗时:"<<(end_time-start_time)/65536.0;
    }
    {
        auto start_time=clock();
        for(unsigned short i(0);i!=end;++i)
        {
            vec x(x0);
            vector<vec> u;
            Ex3::Gauss_Seidel_mat(A,b,x,u);
        }
        auto end_time=clock();
        qDebug()<<"Gauss_Seidel迭代法(矩阵法)耗时:"<<(end_time-start_time)/65536.0;
    }
    {
        auto start_time=clock();
        for(unsigned short i(0);i!=end;++i)
        {
            vec x(x0);
            vector<vec> u;
            Ex3::SOR_ele(A,b,1.08,x,u);
        }
        auto end_time=clock();
        qDebug()<<"SOR法(分量法)耗时:"<<(end_time-start_time)/65536.0;
    }
    {
        auto start_time=clock();
        for(unsigned short i(0);i!=end;++i)
        {
            vec x(x0);
            vector<vec> u;
            Ex3::SOR_mat(A,b,1.08,x,u);
        }
        auto end_time=clock();
        qDebug()<<"SOR法(矩阵法)耗时:"<<(end_time-start_time)/65536.0;
    }
#else
#ifdef DEBUG
    const mat A{{8,1,2},{8,7,2},{4,9,9}},B{{8,10,6,11},{10,16,9,15},{6,9,20,13},{11,15,13,20}};
    const vec b{10,18,17},x0(zeros(3)),b1{1,8,7,9},x1(zeros(4));
    vec x(x1);
    vector<vec> u;
    try {
        Ex3::Jacobi_ele(A,b,x,u);
        Ex3::Jacobi_mat(A,b,x,u);
        Ex3::Gauss_Seidel_ele(A,b,x,u);
        Ex3::Gauss_Seidel_mat(A,b,x,u);
        Ex3::SOR_mat(A,b,1.08,x,u);
        Ex3::SOR_ele(A,b,1.08,x,u);
        Ex3::steepest_descent_method(B,b1,x,u);
        Ex3::conjugate_gradient_method(B,b1,x,u);
        mat H(80,80);
        for(unsigned i(0);i!=80;++i)
            for(unsigned j(0);j!=80;++j)
                H.at(i,j)=1.0/(i+j+1);
        vec e(ones(80)),x(zeros(80));
        vector<mat> t;
        Ex3::Col_Gauss(H,H*e,t,u,1e-6);
        qDebug()<<norm(*(u.end()-1)-e);
        Ex3::SVD_method(H,H*e,x,u);
        qDebug()<<norm(x-e);
    } catch (const char* s) {
        qDebug()<<s;
    } catch (runtime_error&) {
        qDebug()<<"对角线元素为零！";
    }
    for(auto& i:u)_MYFUNCTION::debug_print_vec(i);
#else
#ifdef DEBUG331
    mat A{{7,1,2},{1,8,2},{2,2,9}};
    vec b{10,8,6},x0(zeros(3));
//    Ex3::Jacobi(A,b,x0,"D:/Jacobi.csv");
    Ex3::Gauss_Seidel(A,b,x0,"D:/Gauss_Seidel.csv");
#else
    ui->widget->setTitle("系数矩阵");
    ui->widget_2->setTitle("常数向量");
    ui->widget_3->setTitle("初始向量");
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        double e(calStr(ui->plainTextEdit->toPlainText().toStdString()));
        if(isnan(e)||e<=0)
        {
            QMessageBox::critical(parent,"迭代法","计算精度输入有误！");
            return;
        }
        if(ui->widget->x.empty())
        {
            QMessageBox::critical(parent,"迭代法","系数矩阵为空！");
            return;
        }
        if(ui->widget_2->x.empty())
        {
            QMessageBox::critical(parent,"迭代法","常数向量为空！");
            return;
        }
        if(ui->widget_3->x.empty())
        {
            QMessageBox::critical(parent,"迭代法","初始向量为空！");
            return;
        }
        if(ui->widget->x.n_cols!=ui->widget->x.n_rows)
        {
            QMessageBox::critical(parent,"迭代法","系数矩阵非方阵！");
            return;
        }
        if(ui->widget_2->x.n_rows!=ui->widget->x.n_cols)
        {
            QMessageBox::critical(parent,"迭代法","系数矩阵与常数向量维数不一致！");
            return;
        }
        if(ui->widget_3->x.n_rows!=ui->widget->x.n_cols)
        {
            QMessageBox::critical(parent,"迭代法","系数矩阵与初始向量维数不一致！");
            return;
        }
        if(ui->radioButton->isChecked())
            try {
            Ex3::Jacobi(ui->widget->x,ui->widget_2->x,ui->widget_3->x,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"),e);
        } catch (const char* s) {
            QMessageBox::critical(parent,"迭代法",s);
            return;
        }
        else if(ui->radioButton_2->isChecked())
            try {
            Ex3::Gauss_Seidel(ui->widget->x,ui->widget_2->x,ui->widget_3->x,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"),e);
        } catch (const char* s) {
            QMessageBox::critical(parent,"迭代法",s);
            return;
        }
        else
        try {
            double w(calStr(QInputDialog::getText(parent,"迭代法","请输入松弛因子：").toStdString()));
            if(isnan(w))throw "松弛因子输入有误！";
            Ex3::SOR(ui->widget->x,ui->widget_2->x,w,ui->widget_3->x,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"),e);
        } catch (const char* s) {
            QMessageBox::critical(parent,"迭代法",s);
            return;
        }
        QMessageBox::information(parent,"迭代法","导出成功！");
    });
#endif
#endif
#endif
}

E3t3::~E3t3()
{
    delete ui;
}
