#include "e6t2.h"
#include "ui_e6t2.h"
#include <functional>
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include <QMessageBox>
#include <QInputDialog>
#include <QFileDialog>
#include <QStandardPaths>
#include <QString>
#include <QStringList>
#ifdef DEBUG62
#include <QDebug>
#endif
using namespace arma;
using std::string;
typedef std::function<vec(const double&,const vec&)> func;
typedef std::function<void(const vec&,vec&)> func1;//这里要求结果向量已分配好内存
typedef std::vector<double> vecD;
typedef std::vector<string> vecS;
static double *vec_begin;
static unsigned n;
static FILE* file;
#define EX6T2_THROW_ERROR {fclose(file);throw "区间内存在奇点！";}
#define EX6T2_JUDGE_NAN(x) if(isnan(x))EX6T2_THROW_ERROR
extern double calStr(string);
extern vec calStr_F(const QStringList& S,const vec& x,const vecS& name);
#ifdef DEBUG62
//格式化输出所需的字符串
static char* print_format(nullptr);
//static double** print_content(nullptr);
#ifdef ZSC_ERROR
static char* print_format0(nullptr);
#endif
#endif
namespace _MYFUNCTION {
#ifdef DEBUG62
//设置格式化输出
void set_format(unsigned n)
{
    if(print_format)delete[] print_format;
#ifdef ZSC_ERROR
    if(print_format0)delete[] print_format0;
#endif
    if(!n)
    {
        print_format=nullptr;
#ifdef ZSC_ERROR
        print_format0=nullptr;
#endif
        return;
    }
#ifdef ZSC_ERROR
    unsigned n2(n<<1);
    *(print_format0=new char[n2+1]+n2)='\0';
    *--print_format0='\n';
#endif
    unsigned n1(n*6);
    *(print_format=new char[n1+1]+n1)='\0';
    *--print_format='\n';
    while(--n)
    {
        *--print_format='f';
        *--print_format='4';
        *--print_format='1';
        *--print_format='.';
        *--print_format='%';
        *--print_format=',';
#ifdef ZSC_ERROR
        *--print_format0='0';
        *--print_format0=',';
#endif
    }
    *--print_format='f';
    *--print_format='4';
    *--print_format='1';
    *--print_format='.';
    *--print_format='%';
#ifdef ZSC_ERROR
    *--print_format0='0';
#endif
}
#endif
void print_vec()
{
    double* p(vec_begin);
    fprintf(file,"%.14f",*p);
    unsigned m(n);
    while(--m)fprintf(file,",%.14f",*++p);
}
}
bool isnan(const vec& x)
{
    for(const auto& i:x)
        if(isnan(i))
            return true;
    return false;
}
namespace Ex6 {
/*
 * 显式欧拉法
 * P:文件保存路径
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 * x0:初值点
 *
 * 不检查函数f的维数问题
 */
void Euler(const char* P,const func& f,const double& a,const double& b,const vec& y0,const double& h,double x0)
{
    if(a>=b)throw "区间左端点小于区间右端点！";
    if(x0<a||x0>b)throw "初值点不在区间内！";
    if(!h)throw "步长为零！";
    vec y(y0);
    n=y0.n_elem;
    vec_begin=&y.at(0);
    if(fopen_s(&file,P,"w"))throw "文件保存失败！";
    fprintf(file,"%.14f,",x0);
    _MYFUNCTION::print_vec();
    double x(x0);
    if(h>0)
        while((x0+=h)<=b)
        {
            EX6T2_JUDGE_NAN(y+=h*f(x,y))
            fprintf(file,"\n%.14f,",x=x0);
            _MYFUNCTION::print_vec();
        }
    else
        while((x0+=h)>=a)
        {
            EX6T2_JUDGE_NAN(y+=h*f(x,y))
            fprintf(file,"\n%.14f,",x=x0);
            _MYFUNCTION::print_vec();
        }
    fclose(file);
}
/*
 * 显式欧拉法
 * P:文件保存路径
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 *
 * 不检查函数f的维数问题
 */
inline void Euler(const char* P,const func& f,const double& a,const double& b,const vec& y0,const double& h=0.0078125)
{
    return Euler(P,f,a,b,y0,h,a);
}
/*
 * 4阶Runge-Kutta法
 * P:文件保存路径
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 * x0:初值点
 *
 * 不检查函数f的维数问题
 */
void Runge_Kutta4(const char* P,const func& f,const double& a,const double& b,const vec& y0,const double& h,double x0)
{
    if(a>=b)throw "区间左端点小于区间右端点！";
    if(x0<a||x0>b)throw "初值点不在区间内！";
    if(!h)throw "步长为零！";
    vec y(y0);
    n=y0.n_elem;
    vec_begin=&y.at(0);
    if(fopen_s(&file,P,"w"))throw "文件保存失败！";
    fprintf(file,"%.14f,",x0);
    _MYFUNCTION::print_vec();
    double x(x0),h2(h/2),h6(h2/3);
    if(h>0)
        while((x0+=h)<=b)
        {
            vec k1(f(x,y));
            EX6T2_JUDGE_NAN(k1)
            vec k2(f(x+=h2,y+h2*k1));
            EX6T2_JUDGE_NAN(k2)
            vec k3(f(x,y+h2*k2));
            EX6T2_JUDGE_NAN(k3)
            vec k4(f(x0,y+h*k3));
            EX6T2_JUDGE_NAN(k4)
            y+=h6*(k1+k4+2*(k2+k3));
            fprintf(file,"\n%.14f,",x=x0);
            _MYFUNCTION::print_vec();
        }
    else
        while((x0+=h)>=a)
        {
            vec k1(f(x,y));
            EX6T2_JUDGE_NAN(k1)
            vec k2(f(x+=h2,y+h2*k1));
            EX6T2_JUDGE_NAN(k2)
            vec k3(f(x,y+h2*k2));
            EX6T2_JUDGE_NAN(k3)
            vec k4(f(x0,y+h*k3));
            EX6T2_JUDGE_NAN(k4)
            y+=h6*(k1+k4+2*(k2+k3));
            fprintf(file,"\n%.14f,",x=x0);
            _MYFUNCTION::print_vec();
        }
    fclose(file);
}
/*
 * 4阶Runge-Kutta法
 * P:文件保存路径
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 *
 * 不检查函数f的维数问题
 */
inline void Runge_Kutta4(const char* P,const func& f,const double& a,const double& b,const vec& y0,const double& h=0.0078125)
{
    return Runge_Kutta4(P,f,a,b,y0,h,a);
}
}
E6t2::E6t2(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E6t2)
{
    ui->setupUi(this);
#ifdef DEBUG62
    func f=[](const double&,const vec& x)->vec{
        return vec{4*x.at(0)-2*x.at(0)*x.at(1),x.at(0)*x.at(1)-3*x.at(1)};
    };
    vec y0{2.,3.};
    try {
        Ex6::Euler("D:/Euler.csv",f,0.,10.,y0);
        Ex6::Runge_Kutta4("D:/Runge-Kutta4.csv",f,0.,10.,y0);
    } catch (const char* s) {
        qDebug()<<s;
    }
#else
    ui->widget->setTitle("初值");
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        try {
            double a(calStr(ui->lineEdit->text().toStdString()));
            if(isnan(a))throw "区间左端点输入有误！";
            double b(calStr(ui->lineEdit_2->text().toStdString()));
            if(isnan(b))throw "区间右端点输入有误！";
            double h(calStr(ui->lineEdit_3->text().toStdString()));
            if(isnan(h))throw "步长输入有误！";
            vec& x0(ui->widget->x);
            if(x0.empty())throw "初始向量为空！";
            QStringList S(ui->plainTextEdit->toPlainText().split('\n'));
            if(int(x0.n_elem)!=S.size())throw "函数维数与初值向量维数不一致！";
            double t0(calStr(QInputDialog::getText(this,"常微分方程组数值解法","请输入初值点：").toStdString()));
            if(isnan(t0))throw "初值点输入有误！";
            vecS name(x0.n_elem);
            auto p(name.begin());
            for(unsigned i(0);i!=x0.n_elem;)
                *p++=(QInputDialog::getText(this,"常微分方程组数值解法","请输入第"+QString::number(++i)+"个变量名：").toStdString());
            name.push_back("t");
            func f=[=](const double& t,const vec& x)->vec{
                vec _x(x);
                _x.resize(x.n_elem+1).at(x.n_elem)=t;
                return calStr_F(S,_x,name);
            };
            const char* P(QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()").toLatin1().data());
            if(ui->radioButton->isChecked())
                Ex6::Euler(P,f,a,b,x0,h,t0);
            else
                Ex6::Runge_Kutta4(P,f,a,b,x0,h,t0);
            QMessageBox::information(this,"常微分方程组数值解法","导出成功！");
        } catch (const char* s) {
            QMessageBox::critical(this,"常微分方程组数值解法",s);
        }
    });
#endif
}

E6t2::~E6t2()
{
    delete ui;
#ifdef DEBUG62
    if(print_format)delete[] print_format;
#ifdef ZSC_ERROR
    if(print_format0)delete[] print_format0;
#endif
#endif
}
