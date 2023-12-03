#include "e2t8.h"
#include "ui_e2t8.h"
#include <functional>
#include <QFile>
#include <QStandardPaths>
#include <QFileDialog>
#include <math.h>
#include <vector>
#include <string>
#include <QMessageBox>
#ifdef DEBUG28
#include <QDebug>
#endif
using namespace std;
extern double calStr(string);
extern double calStr_x(string,const double&);
namespace _MYFUNCTION {
//符号函数
template <class T>
inline int sgn(const T& x)noexcept
{
    if(x==0)return 0;
    if(x>0)return 1;
    return -1;
}
}
#ifdef DEBUG28
namespace Ex1 {
//template <class T>
//T &iterative_method(T& x,std::function<void(T&)> func,std::function<bool(const T&)> exit,std::function<void(const T&)> out);
/*
 * 迭代法(模板)
 * x   : 迭代自变量
 * func: 迭代函数
 * exit: 终止标准
 *  true  - 继续迭代
 *  false - 输出结果
 *  第1个参数: x_n
 *  第2个参数: x_{n-1}
 */
template <class T>
T &iterative_method(T& x,const function<void(T&)>& func,const function<bool(T,T)>& exit,const function<void(const double&)>& out)
{
    T t(x);
    out(x);
    func(x);
    out(x);
    while(exit(x,t))
    {
        t=x;
        func(x);
        out(x);
    }
    return x;
}
}
#endif
#if defined(DEBUG28)&&!defined(DEBUG282)
unsigned secant_num(0);
#endif
namespace Ex2 {
/*
 * 割线法
 * f   : 原方程函数
 * x1  : 第1个值
 * x2  : 第2个值
 * e   : 精度
 * path: 迭代保存路径
 */
void secant_method(const function<double(double)>& f,double x1,double x2,double e=1e-6,const QString& path=QStandardPaths::writableLocation(QStandardPaths::DesktopLocation))
{
    QFile file(path);
    if(!(file.open(QIODevice::WriteOnly)))throw "文件保存失败！";
    double f1(f(x1)),f2;
    if(isnan(f1))throw "区间内存在奇点！";
    file.write((to_string(x1)+'\n'+to_string(x2)).c_str());
    while(abs(x1-x2)>=e)
    {
        if(isnan(f2=f(x2)))
        {
            file.close();
            throw "区间内存在奇点！";
        }
        double t=x2-(x2-x1)*f2/(f2-f1);
        file.write(('\n'+to_string(t)).c_str());
        x1=x2;
        x2=t;
        f1=f2;
#if defined(DEBUG28)&&!defined(DEBUG282)
        ++secant_num;
#endif
    }
#if defined(DEBUG28)&&!defined(DEBUG282)
    qDebug()<<f(x2);
#endif
    file.close();
}
/*
 * 抛物线法
 * f   : 原方程函数
 * x1  : 第1个值
 * x2  : 第2个值
 * x3  : 第3个值
 * e   : 精度
 * path: 迭代保存路径
 */
void parabolic_method(const function<double(double)>& f,double x1,double x2,double x3,double e=1e-6,const QString& path=QStandardPaths::writableLocation(QStandardPaths::DesktopLocation))
{
    QFile file(path);
    if(!(file.open(QIODevice::WriteOnly)))throw "文件保存失败！";
    double f1(f(x1)),f2(f(x2)),f3,omega0((f2-f1)/(x2-x1));
    if(isnan(f1)||isnan(f2))
    {
        file.close();
        throw "区间内存在奇点！";
    }
    file.write((to_string(x1)+'\n'+to_string(x2)+'\n'+to_string(x3)).c_str());
    while(abs(x1-x2)>=e)
    {
        if(isnan(f3=f(x3)))
        {
            file.close();
            throw "区间内存在奇点！";
        }
        double omega1=(f3-f2)/(x3-x2),d123=(omega1-omega0)/(x3-x1),omega=omega1+(x3-x2)*d123,delta=omega*omega-4*f3*d123;
        if(delta<0)
        {
            file.close();
            throw "迭代出现复根！";
        }
        double t=x3-2*f3/(omega+_MYFUNCTION::sgn(omega)*sqrt(delta));
        file.write(('\n'+to_string(t)).c_str());
        x1=x2;
        x2=x3;
        x3=t;
        f1=f2;
        f2=f3;
        omega0=omega1;
    }
    file.close();
}
/*
 * 逆二次插值法
 * f   : 原方程函数
 * x1  : 第1个值
 * x2  : 第2个值
 * x3  : 第3个值
 * e   : 精度
 * path: 迭代保存路径
 */
void inverse_quadratic_interpolation(const function<double(double)>& f,double x1,double x2,double x3,double e=1e-6,const QString& path=QStandardPaths::writableLocation(QStandardPaths::DesktopLocation))
{
    QFile file(path);
    if(!(file.open(QIODevice::WriteOnly)))throw "文件保存失败！";
    double f1(f(x1)),f2(f(x2)),f3,omega0((x2-x1)/(f2-f1));
    if(isnan(f1)||isnan(f2))
    {
        file.close();
        throw "区间内存在奇点！";
    }
    file.write((to_string(x1)+'\n'+to_string(x2)+'\n'+to_string(x3)).c_str());
    while(abs(x1-x2)>=e)
    {
        if(isnan(f3=f(x3)))
        {
            file.close();
            throw "区间内存在奇点！";
        }
        double omega1=(x3-x2)/(f3-f2),t=x3-omega1*f3+(omega1-omega0)*f3*f2/(f3-f1);
        file.write(('\n'+to_string(t)).c_str());
        x1=x2;
        x2=x3;
        x3=t;
        f1=f2;
        f2=f3;
        omega0=omega1;
    }
    file.close();
}
/*
 * 布伦特法
 * f   : 原方程函数
 * a   : 区间左端点
 * b   : 区间右端点
 * e   : 精度
 * path: 迭代保存路径
 */
void Brent_method(const function<double(double)>& f,double a,double b,double e=1e-6,const QString& path=QStandardPaths::writableLocation(QStandardPaths::DesktopLocation))
{
    QFile file(path);
    if(!(file.open(QIODevice::WriteOnly)))throw "文件保存失败！";
    double fa(f(a)),fb(f(b));
    if(isnan(fa)||isnan(fb))
    {
        file.close();
        throw "区间内存在奇点！";
    }
    if(abs(fa)<e)
    {
        file.write(to_string(fa).c_str());
        file.close();
        return;
    }
    if(abs(fb)<e)
    {
        file.write(to_string(fb).c_str());
        file.close();
        return;
    }
    bool faG0(fa>0),fbG0(fb>0);
    if(faG0&&fbG0||!faG0&&!fbG0)
    {
        file.close();
        return parabolic_method(f,a,(a+b)/2,b,e,path);
    }
    file.write((to_string(a)+'\n'+to_string(b)).c_str());
    double pre_m1(a),pre_fm1(fa),pre_m2(b),pre_fm2(fb),pre_dx(b-a),pre_df(fb-fa),pre_dfx(pre_dx/pre_df);
    do
    {
        double m((a+b)/2),fm(f(m)),dx(m-pre_m2),df(fm-pre_fm2),x1,dfx(dx/df);
        if(isnan(fm)||abs(fm)<e||abs(fm-fa)<e||abs(fm-fb)<e)
        {
            file.close();
            return;
        }
        if(faG0)
            if(fm>0)
            {
                a=m;
                fa=fm;
            }
            else
            {
                b=m;
                fb=fm;
            }
        else
            if(fm>0)
            {
                b=m;
                fb=fm;
            }
            else
            {
                a=m;
                fa=fm;
            }
        if(!df||!pre_df)
        {
            pre_m1=a;
            pre_fm1=fa;
            pre_m2=b;
            pre_fm2=fb;
            m=(a+b)/2;
            fm=f(m);
        }
        if(pre_fm1!=fm)
            x1=m-fm*dfx+pre_m2*fm*(dfx-pre_dfx)/(fm-pre_fm1);
        else
            if(df)
                x1=m-fm*dfx;
            else
                x1=pre_m2-pre_fm2*pre_dfx;
        double fx1(f(x1));
        if(isnan(fx1))
        {
            file.close();
            throw "区间内存在奇点！";
        }
        if((fx1-fa<0)^(fx1-fb<0))x1=(a+b)/2;
        pre_dfx=(pre_dx=x1-(pre_m1=pre_m2))/(pre_df=fx1-(pre_fm1=pre_fm2));
        pre_fm2=fx1;
        file.write(('\n'+to_string(pre_m2=x1)).c_str());
    }while(true);
}
}
E2t8::E2t8(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E2t8)
{
    ui->setupUi(this);
#ifdef DEBUG28
#ifdef DEBUG282
    auto f=[](double x)->double{
        return ((x+2)*x+10)*x-20;
    };
    auto Newton_f=[](double& x){
        x-=(((x+2)*x+10)*x-20)/((3*x+4)*x+10);
    };//Newton迭代函数, 使用ggb推导, 请参考Newton iteration.ggb
    auto Newton_exit=[](const double& x,const double& y)->bool{
        return abs(x-y)>=1e-6;
    };
    FILE* file;
    try {
#ifdef DEBUG281
        int8_t i(-128),end(i);
        char s[14];
        FILE* total_file;
        if(fopen_s(&total_file,"D:/1/total.csv","w"))throw -129;
        do
        {
            double x0(i);
            sprintf(s,"D:/1/%hhd.csv",i);
            if(fopen_s(&file,s,"w"))throw i;
            unsigned k(0);
            auto Newton_out=[file,&k](const double& x){
                fprintf(file,"%u,%.14f\n",++k,x);
            };
            double error(f(Ex1::iterative_method<double>(x0,Newton_f,Newton_exit,Newton_out)));
            fprintf(total_file,"%hhd,%.14f,%u\n",i,error,k);
            fclose(file);
        }while(++i!=end);
        fclose(total_file);
#else
        int8_t i(-128),end(i);
        char s[20];
        FILE* total_file;
        if(fopen_s(&total_file,"D:/2/total.csv","w"))throw -129;
        do
        {
            double x0(1.37+i/128.),x(x0);
            sprintf(s,"D:/2/%g.csv",x0);
            if(fopen_s(&file,s,"w"))throw i;
            unsigned k(0);
            auto Newton_out=[file,&k](const double& x){
                fprintf(file,"%u,%.14f\n",++k,x);
            };
            double error(f(Ex1::iterative_method<double>(x,Newton_f,Newton_exit,Newton_out)));
            fprintf(total_file,"%f,%.14f,%u\n",x0,error,k);
            fclose(file);
        }while(++i!=end);
        fclose(total_file);
#endif
    } catch (uint8_t i) {
        qDebug()<<(int)i;
    } catch (int i) {
        qDebug()<<i;
    }
#else
    auto f=[](double x)->double{
        return (x*x-2)*x-5;
    };
    auto Newton_f=[](double& x){
        x-=((x*x-2)*x-5)/(3*x*x-2);
    };
    auto Newton_exit=[](const double& x,const double& y)->bool{
        return abs(x-y)>=1e-6;
    };
    FILE* file;
    try {
        double x0(4),x1(3.8);
        Ex2::secant_method(f,x0,x1,1e-6,"D:/secant.csv");
        qDebug()<<secant_num;
        if(fopen_s(&file,"D:/Newton.csv","w"))throw 0;
        unsigned Newton_num(-1);
        auto Newton_out=[file,&Newton_num](const double& x){
            fprintf(file,"%.14f\n",x);
            ++Newton_num;
        };
        qDebug()<<f(Ex1::iterative_method<double>(x1,Newton_f,Newton_exit,Newton_out));
        qDebug()<<Newton_num;
    } catch (int i) {
        qDebug()<<i;
    } catch (const char* s) {
        qDebug()<<s;
    }
#endif
#else
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        try {
            double x1(calStr(ui->lineEdit_2->text().toStdString()));
            if(isnan(x1))throw "x1输入有误！";
            double x2(calStr(ui->lineEdit_3->text().toStdString()));
            if(isnan(x2))throw "x2输入有误！";
            double e(calStr(ui->lineEdit_5->text().toStdString()));
            if(isnan(e))throw "计算精度输入有误！";
            string s(ui->lineEdit->text().toStdString());
            auto f=[s](double x)->double{
                return calStr_x(s,x);
            };
            if(ui->radioButton->isChecked())
                Ex2::secant_method(f,x1,x2,e,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"));
            else if(ui->radioButton_4->isChecked())
                Ex2::Brent_method(f,x1,x2,e,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"));
            else
            {
                double x3(calStr(ui->lineEdit_4->text().toStdString()));
                if(isnan(x3))throw "x3输入有误！";
                if(ui->radioButton_2->isChecked())
                    Ex2::parabolic_method(f,x1,x2,x3,e,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"));
                else
                    Ex2::inverse_quadratic_interpolation(f,x1,x2,x3,e,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"));
            }
            QMessageBox::information(this,"迭代法(选择)","导出成功！");
        } catch (const char* s) {
            QMessageBox::critical(this,"迭代法(选择)",s);
        }
    });
#endif
}

E2t8::~E2t8()
{
    delete ui;
}
