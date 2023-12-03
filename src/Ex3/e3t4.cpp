#include "e3t4.h"
#include "ui_e3t4.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QStandardPaths>
using namespace arma;
#define THROW_DIV_0 throw "方阵为奇异矩阵！";
#define THROW_DIV_0_DEL {free(B);THROW_DIV_0}
//#define THROW_DIV_0_ALL {free(A);free(B);free(C);THROW_DIV_0}
#define JUDGE_0(x) if(!(x))THROW_DIV_0
#define JUDGE_0_DEL(x) if(x)THROW_DIV_0_DEL
//#define JUDGE_0_ALL(x) if(x)THROW_DIV_0_ALL
#ifdef DEBUG34
#include <QDebug>
#endif
namespace _MYFUNCTION {
//追赶过程, 要求b[0]不为零, n>1
//结果输出在d中
bool chase(unsigned n,const double* a,double* b,const double* c,double* d)
{
    unsigned k(n);
    while(--k)
    {
        if(!*b)return true;
        double m(*a++/ *b),&pre_d(*d++);
        *++b-=m**c++;
        *d-=m*pre_d;
    }
    if(!*b)return true;
    *d/=*b;
    while(--n)
    {
        double &pre_d(*d--);
        (*d-=*--c*pre_d)/=*--b;
    }
    return false;
}
//b[0]为零的追赶, n>2
bool chase_0(unsigned n,const double* a,double* b,const double* c,double* d)
{
    if(!*c||!*a)return true;
    const double t(*d/ *c++);
    *++d-=t**(++b)++;
    *++d-=t**++a;
    unsigned k(n);
    while(--k)
    {
        if(!*b)return true;
        const double m(*++a/ *b),&pre_d(*d++);
        *++b-=m**++c;
        *d-=m*pre_d;
    }
    if(!*b)return true;
    *d/=*b;
    a-=n;
    while(--n)
    {
        const double &pre_x(*d--);
        (*d-=*c--*pre_x)/=*--b;
    }
    const double &x3(*d);
    double &x2(*--d);
    *--d=(x2-*c*x3)/ *a;
    x2=t;
    return false;
}
}
namespace Ex3 {
/*
 * 追赶法
 * R:解向量
 * a:左下方的副对角线
 * b:主对角线
 * c:右上方的副对角线
 * d:常数项
 */
void Thomas_algorithm(vec& R,const vec& a,const vec& b,const vec& c,const vec& d)
{
    if(b.n_elem!=d.n_elem)throw "系数矩阵阶数与常数向量维数不相等！";
    if(a.n_elem!=b.n_elem-1||a.n_elem!=c.n_elem)throw "主、副对角线维数不匹配！";
    if(a.empty())
    {
        double t(d.at(0)/b.at(0));
        JUDGE_0(t)
        R={t};
        return;
    }
    unsigned n(b.n_elem);
    R.set_size(n);
    const double *Bp(&b.at(--n)),*Dp(&d.at(n));
    double *B((double*)malloc(sizeof(double)*b.n_elem)),*bp(B+n),*dp(&R.at(n));
    *bp=*Bp,*dp=*Dp;
    do *--bp=*--Bp,*--dp=*--Dp;while(--n);
    if(*bp)//第1个元素不为0, 可以直接进入追赶过程
    {
        JUDGE_0_DEL(_MYFUNCTION::chase(b.n_elem,&a.at(0),B,&c.at(0),dp))
        free(B);
        return;
    }
    if(!--(n=a.n_elem))//系数矩阵只有2阶
    {
        const double &C(c.at(0)),&A(a.at(0));
        double &ele_2(*(dp+1));
        JUDGE_0_DEL(!C||!A)
        *dp=(Dp[1]-B[1]*(ele_2=*Dp/C))/A;
        free(B);
        return;
    }
//    const unsigned A_size(sizeof(double)*n);
//    const double *AP(&a.at(0)),*CP(&c.at(0));
//    JUDGE_0_DEL(!*CP)
//    double *A((double*)malloc(A_size)),*C((double*)malloc(A_size)),last_sol(*Dp/ *CP);
    JUDGE_0_DEL(_MYFUNCTION::chase_0(n,&a.at(0),B,&c.at(0),dp))
//    free(A);
    free(B);
//    free(C);
}
}
E3t4::E3t4(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E3t4)
{
    ui->setupUi(this);
#ifdef DEBUG34
    vec R,a{1,1,1},c(a),b(4,fill::value(-4)),d(ones(4));
    try {
        Ex3::Thomas_algorithm(R,a,b,c,d);
        R.print("res1=");
        b.at(0)=0;
        Ex3::Thomas_algorithm(R,a,b,c,d);
        R.print("res2=");
    } catch (const char* s) {
        qDebug()<<s;
    }
#endif
    ui->widget->setTitle("左下方的副对角线");
    ui->widget_2->setTitle("主对角线");
    ui->widget_3->setTitle("右上方的副对角线");
    ui->widget_4->setTitle("常数项");
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        try {
            vec R;
            Ex3::Thomas_algorithm(R,ui->widget->x,ui->widget_2->x,ui->widget_3->x,ui->widget_4->x);
            if(!R.save(QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()").toStdString()))throw "文件保存失败！";
            QMessageBox::information(this,"追赶法","文件导出成功！");
        } catch (const char* s) {
            QMessageBox::critical(this,"追赶法",s);
        }
    });
}

E3t4::~E3t4()
{
    delete ui;
}
