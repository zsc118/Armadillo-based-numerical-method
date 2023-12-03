#include "e6t1.h"
#include "ui_e6t1.h"
#include <functional>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <string>
#include <QMessageBox>
#include <QFileDialog>
#include <QStandardPaths>
typedef std::function<double(const double&,const double&)> func;
typedef std::function<double(const double&)> func1;
typedef std::vector<double> vec;
using std::string;
extern double calStr_xy(string,const double&,const double&);
extern double calStr(string);
#if defined(ZSC_ERROR)&&defined(DEBUG61)
func1 ans=[](const double& x)->double{return exp(-2*x)+x*x;};
#endif
#ifdef DEBUG61
#include <QDebug>
#endif
#ifdef DEBUG611
#include <QDebug>
func1 ans=[](const double& x)->double{return exp(-50*x)/2;};
#endif
namespace Ex6 {
/*
 * 显式欧拉法
 * R:结果向量
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 * x0:初值点
 */
void Euler(vec& R,const func& f,const double& a,const double& b,const double& y0,const double& h,double x0)
{
    if(a>=b)throw "区间左端点小于区间右端点！";
    if(x0<a||x0>b)throw "初值点不在区间内！";
    if(!h)throw "步长为零！";
    R={y0};
    double y(y0),x(x0);
    if(h>0)
        while((x0+=h)<=b)
        {
            if(isnan(y+=h*f(x,y)))
                throw "区间内存在奇点！";
            R.push_back(y);
            x=x0;
        }
    else
        while((x0+=h)>=a)
        {
            if(isnan(y+=h*f(x,y)))
                throw "区间内存在奇点！";
            R.push_back(y);
            x=x0;
        }
}
/*
 * 显式欧拉法
 * R:结果向量
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 */
inline void Euler(vec& R,const func& f,const double& a,const double& b,const double& y0,const double& h=0.0078125)
{
    return Euler(R,f,a,b,y0,h,a);
}
/*
 * 显式欧拉法
 * P:结果保存路径
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 * x0:初值点
 */
void Euler(const char* P,const func& f,const double& a,const double& b,const double& y0,const double& h,double x0)
{
    if(a>=b)throw "区间左端点小于区间右端点！";
    if(x0<a||x0>b)throw "初值点不在区间内！";
    if(!h)throw "步长为零！";
    FILE* file;
    if(fopen_s(&file,P,"w"))throw "文件保存失败！";
#ifdef ZSC_ERROR
    fprintf(file,"0,0\n");
#else
    fprintf(file,"%.14f,%.14f\n",x0,y0);
#endif
    double y(y0),x(x0);
    if(h>0)
        while((x0+=h)<=b)
        {
            if(isnan(y+=h*f(x,y)))
            {
                fclose(file);
                throw "区间内存在奇点！";
            }
#ifdef ZSC_ERROR
            fprintf(file,"%.14f,%.14f\n",x=x0,y-ans(x0));
#else
            fprintf(file,"%.14f,%.14f\n",x=x0,y);
#endif
        }
    else
        while((x0+=h)>=a)
        {
            if(isnan(y+=h*f(x,y)))
            {
                fclose(file);
                throw "区间内存在奇点！";
            }
#ifdef ZSC_ERROR
            fprintf(file,"%.14f,%.14f\n",x=x0,y-ans(x0));
#else
            fprintf(file,"%.14f,%.14f\n",x=x0,y);
#endif
        }
    fclose(file);
}
/*
 * 显式欧拉法
 * P:结果保存路径
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 * x0:初值点
 */
inline void Euler(const char* P,const func& f,const double& a,const double& b,const double& y0,const double& h=0.0078125)
{
    return Euler(P,f,a,b,y0,h,a);
}
/*
 * 隐式欧拉法
 * R:结果向量
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 * x0:初值点
 */
void imEuler(vec& R,const func& f,const double& a,const double& b,const double& y0,const double& h,double x0)
{
    if(a>=b)throw "区间左端点小于区间右端点！";
    if(x0<a||x0>b)throw "初值点不在区间内！";
    if(!h)throw "步长为零！";
    R={y0};
    double y(y0),x(x0);
    if(h>0)
        while((x0+=h)<=b)
        {
            double t(y+h*f(x,y));
            if(isnan(t)||isnan(y+=h*f(x=x0,t)))
                throw "区间内存在奇点！";
            R.push_back(y);
        }
    else
        while((x0+=h)>=a)
        {
            double t(y+h*f(x,y));
            if(isnan(t)||isnan(y+=h*f(x=x0,t)))
                throw "区间内存在奇点！";
            R.push_back(y);
        }
}
/*
 * 隐式欧拉法
 * R:结果向量
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 */
inline void imEuler(vec& R,const func& f,const double& a,const double& b,const double& y0,const double& h=0.0078125)
{
    return imEuler(R,f,a,b,y0,h,a);
}
/*
 * 隐式欧拉法
 * P:结果保存路径
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 * x0:初值点
 */
void imEuler(const char* P,const func& f,const double& a,const double& b,const double& y0,const double& h,double x0)
{
    if(a>=b)throw "区间左端点小于区间右端点！";
    if(x0<a||x0>b)throw "初值点不在区间内！";
    if(!h)throw "步长为零！";
    FILE* file;
    if(fopen_s(&file,P,"w"))throw "文件保存失败！";
#ifdef ZSC_ERROR
    fprintf(file,"0,0\n");
#else
    fprintf(file,"%.14f,%.14f\n",x0,y0);
#endif
    double y(y0),x(x0);
    if(h>0)
        while((x0+=h)<=b)
        {
            double t(y+h*f(x,y));
            if(isnan(t)||isnan(y+=h*f(x0,t)))
            {
                fclose(file);
                throw "区间内存在奇点！";
            }
#ifdef ZSC_ERROR
            fprintf(file,"%.14f,%.14f\n",x=x0,y-ans(x0));
#else
            fprintf(file,"%.14f,%.14f\n",x=x0,y);
#endif
        }
    else
        while((x0+=h)>=a)
        {
            double t(y+h*f(x,y));
            if(isnan(t)||isnan(y+=h*f(x0,t)))
            {
                fclose(file);
                throw "区间内存在奇点！";
            }
#ifdef ZSC_ERROR
            fprintf(file,"%.14f,%.14f\n",x=x0,y-ans(x0));
#else
            fprintf(file,"%.14f,%.14f\n",x=x0,y);
#endif
        }
    fclose(file);
}
/*
 * 隐式欧拉法
 * P:结果保存路径
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 * x0:初值点
 */
inline void imEuler(const char* P,const func& f,const double& a,const double& b,const double& y0,const double& h=0.0078125)
{
    return imEuler(P,f,a,b,y0,h,a);
}
/*
 * 4阶Runge-Kutta法
 * P:结果保存路径
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 * x0:初值点
 */
void Runge_Kutta4(const char* P,const func& f,const double& a,const double& b,const double& y0,const double& h,double x0)
{
    if(a>=b)throw "区间左端点小于区间右端点！";
    if(x0<a||x0>b)throw "初值点不在区间内！";
    if(!h)throw "步长为零！";
    FILE* file;
    if(fopen_s(&file,P,"w"))throw "文件保存失败！";
#ifdef ZSC_ERROR
    fprintf(file,"0,0\n");
#else
    fprintf(file,"%.14f,%.14f\n",x0,y0);
#endif
    double y(y0),x(x0);
    const double h2(h/2),h6(h2/3);
    if(h>0)
        while((x0+=h)<=b)
        {
            double k1(f(x,y));
            if(isnan(k1))
            {
                fclose(file);
                throw "区间内存在奇点！";
            }
            double k2(f(x+=h2,y+k1*h2));
            if(isnan(k2))
            {
                fclose(file);
                throw "区间内存在奇点！";
            }
            double k3(f(x,y+k2*h2));
            if(isnan(k3))
            {
                fclose(file);
                throw "区间内存在奇点！";
            }
            double k4(f(x0,y+h*k3));
            if(isnan(k4))
            {
                fclose(file);
                throw "区间内存在奇点！";
            }
#ifdef ZSC_ERROR
            fprintf(file,"%.14f,%.14f\n",x=x0,(y+=h6*(k1+k4+2*(k2+k3)))-ans(x0));
#else
            fprintf(file,"%.14f,%.14f\n",x=x0,y+=h6*(k1+k4+2*(k2+k3)));
#endif
        }
    else
        while((x0+=h)>=a)
        {
            double k1(f(x,y));
            if(isnan(k1))
            {
                fclose(file);
                throw "区间内存在奇点！";
            }
            double k2(f(x+=h2,y+k1*h2));
            if(isnan(k2))
            {
                fclose(file);
                throw "区间内存在奇点！";
            }
            double k3(f(x,y+k2*h2));
            if(isnan(k3))
            {
                fclose(file);
                throw "区间内存在奇点！";
            }
            double k4(f(x0,y+h*k3));
            if(isnan(k4))
            {
                fclose(file);
                throw "区间内存在奇点！";
            }
#ifdef ZSC_ERROR
            fprintf(file,"%.14f,%.14f\n",x=x0,(y+=h6*(k1+k4+2*(k2+k3)))-ans(x0));
#else
            fprintf(file,"%.14f,%.14f\n",x=x0,y+=h6*(k1+k4+2*(k2+k3)));
#endif
        }
    fclose(file);
}
/*
 * 4阶Runge-Kutta法
 * P:结果保存路径
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 * x0:初值点
 */
inline void Runge_Kutta4(const char* P,const func& f,const double& a,const double& b,const double& y0,const double& h=0.0078125)
{
    return Runge_Kutta4(P,f,a,b,y0,h,a);
}
/*
 * 4阶Runge-Kutta法
 * R:结果向量
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 * x0:初值点
 */
void Runge_Kutta4(vec& R,const func& f,const double& a,const double& b,const double& y0,const double& h,double x0)
{
    if(a>=b)throw "区间左端点小于区间右端点！";
    if(x0<a||x0>b)throw "初值点不在区间内！";
    if(!h)throw "步长为零！";
    R={y0};
    double y(y0),x(x0);
    const double h2(h/2),h6(h2/3);
    if(h>0)
        while((x0+=h)<=b)
        {
            double k1(f(x,y));
            if(isnan(k1))
                throw "区间内存在奇点！";
            double k2(f(x+=h2,y+k1*h2));
            if(isnan(k2))
                throw "区间内存在奇点！";
            double k3(f(x,y+k2*h2));
            if(isnan(k3))
                throw "区间内存在奇点！";
            double k4(f(x0,y+h*k3));
            if(isnan(k4))
                throw "区间内存在奇点！";
            x=x0;
            R.push_back(y+=h6*(k1+k4+2*(k2+k3)));
        }
    else
        while((x0+=h)>=a)
        {
            double k1(f(x,y));
            if(isnan(k1))
                throw "区间内存在奇点！";
            double k2(f(x+=h2,y+k1*h2));
            if(isnan(k2))
                throw "区间内存在奇点！";
            double k3(f(x,y+k2*h2));
            if(isnan(k3))
                throw "区间内存在奇点！";
            double k4(f(x0,y+h*k3));
            if(isnan(k4))
                throw "区间内存在奇点！";
            x=x0;
            R.push_back(y+=h6*(k1+k4+2*(k2+k3)));
        }
}
/*
 * 4阶Runge-Kutta法
 * R:结果向量
 * f:显式微分方程的右端函数
 * a:区间左端点
 * b:区间右端点
 * y0:解在x0处的初值
 * h:步长
 */
inline void Runge_Kutta4(vec& R,const func& f,const double& a,const double& b,const double& y0,const double& h=0.0078125)
{
    return Runge_Kutta4(R,f,a,b,y0,h,a);
}
}
#if defined(DEBUG61)||defined(DEBUG611)
#ifndef ZSC_ERROR
void print_ans(const char* p,const func1& f,double a,const double& b,unsigned n=10000)
{
    FILE* file;
    if(fopen_s(&file,p,"w"))throw "答案保存失败！";
    double h((b-a)/n);
    do
    {
        fprintf(file,"%.14f,%.14f\n",a,f(a));
        a+=h;
    }while(--n);
    fprintf(file,"%.14f,%.14f",a,f(a));
}
#endif
#endif
E6t1::E6t1(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E6t1)
{
    ui->setupUi(this);
#ifdef DEBUG61
    func f=[](const double& x,const double& y)->double{return 2*(x*(x+1)-y);};
#ifndef ZSC_ERROR
    func1 ans=[](const double& x)->double{return exp(-2*x)+x*x;};
#endif
    try {
#ifdef ZSC_ERROR
        Ex6::Euler("D:/0.1err.csv",f,0.,1.,1.,0.1);
        Ex6::Euler("D:/0.05err.csv",f,0.,1.,1.,0.05);
        Ex6::Euler("D:/0.01err.csv",f,0.,1.,1.,0.01);
        Ex6::imEuler("D:/0.1err.csv",f,0.,1.,1.,0.1);
        Ex6::imEuler("D:/0.05err.csv",f,0.,1.,1.,0.05);
        Ex6::imEuler("D:/0.01err.csv",f,0.,1.,1.,0.01);
        Ex6::Runge_Kutta4("D:/0.1err.csv",f,0.,1.,1.,0.1);
        Ex6::Runge_Kutta4("D:/0.05err.csv",f,0.,1.,1.,0.05);
        Ex6::Runge_Kutta4("D:/0.01err.csv",f,0.,1.,1.,0.01);
#else
        Ex6::Euler("D:/0.1.csv",f,0.,1.,1.,0.1);
        Ex6::Euler("D:/0.05.csv",f,0.,1.,1.,0.05);
        Ex6::Euler("D:/0.01.csv",f,0.,1.,1.,0.01);
        Ex6::imEuler("D:/0.1.csv",f,0.,1.,1.,0.1);
        Ex6::imEuler("D:/0.05.csv",f,0.,1.,1.,0.05);
        Ex6::imEuler("D:/0.01.csv",f,0.,1.,1.,0.01);
        Ex6::Runge_Kutta4("D:/0.1.csv",f,0.,1.,1.,0.1);
        Ex6::Runge_Kutta4("D:/0.05.csv",f,0.,1.,1.,0.05);
        Ex6::Runge_Kutta4("D:/0.01.csv",f,0.,1.,1.,0.01);
        print_ans("D:/ans.csv",ans,0,1.);
#endif
    } catch (const char* s) {
        qDebug()<<s;
    }
#elif defined(DEBUG611)
    func f=[](const double&,const double& y)->double{return -50*y;};
    try {
//        Ex6::Euler("D:/Euler.csv",f,0.,1.,0.5);
//        Ex6::imEuler("D:/imEuler.csv",f,0.,1.,0.5);
//        Ex6::Runge_Kutta4("D:/Runge-Kutta.csv",f,0.,1.,0.5);
//        print_ans("D:/ans.csv",ans,0,1.);
//        Ex6::Euler("D:/8.csv",f,0.,1.,0.5,0.125);
//        Ex6::Euler("D:/16.csv",f,0.,1.,0.5,0.0625);
//        Ex6::Euler("D:/32.csv",f,0.,1.,0.5,0.03125);
//        Ex6::Euler("D:/64.csv",f,0.,1.,0.5,0.015625);
//        Ex6::Euler("D:/128.csv",f,0.,1.,0.5);
//        Ex6::imEuler("D:/8.csv",f,0.,1.,0.5,0.125);
//        Ex6::imEuler("D:/16.csv",f,0.,1.,0.5,0.0625);
//        Ex6::imEuler("D:/32.csv",f,0.,1.,0.5,0.03125);
//        Ex6::imEuler("D:/64.csv",f,0.,1.,0.5,0.015625);
//        Ex6::imEuler("D:/128.csv",f,0.,1.,0.5);
        Ex6::Runge_Kutta4("D:/8.csv",f,0.,1.,0.5,0.125);
        Ex6::Runge_Kutta4("D:/16.csv",f,0.,1.,0.5,0.0625);
        Ex6::Runge_Kutta4("D:/32.csv",f,0.,1.,0.5,0.03125);
        Ex6::Runge_Kutta4("D:/64.csv",f,0.,1.,0.5,0.015625);
        Ex6::Runge_Kutta4("D:/128.csv",f,0.,1.,0.5);
    } catch (const char* s) {
        qDebug()<<s;
    }
#else
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        try {
            double a(calStr(ui->lineEdit_2->text().toStdString()));
            if(isnan(a))throw "自变量范围左端点输入有误！";
            double b(calStr(ui->lineEdit_3->text().toStdString()));
            if(isnan(b))throw "自变量范围右端点输入有误！";
            double x0(calStr(ui->lineEdit_4->text().toStdString()));
            if(isnan(x0))throw "初值点输入有误！";
            double h(calStr(ui->lineEdit_5->text().toStdString()));
            if(isnan(h))throw "步长输入有误！";
            double y0(calStr(ui->lineEdit_6->text().toStdString()));
            if(isnan(y0))throw "初值输入有误！";
            const char* p(QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()").toLatin1().data());
            string s(ui->lineEdit->text().toStdString());
            func f=[s](const double& x,const double& y)->double{
                return calStr_xy(s,x,y);
            };
            if(ui->radioButton->isChecked())
                Ex6::Euler(p,f,a,b,y0,h,x0);
            else if(ui->radioButton_2->isChecked())
                Ex6::imEuler(p,f,a,b,y0,h,x0);
            else
                Ex6::Runge_Kutta4(p,f,a,b,y0,h,x0);
            QMessageBox::information(this,"常微分方程数值解法","导出成功！");
        } catch (const char* s) {
            QMessageBox::critical(this,"常微分方程数值解法",s);
        }
    });
#endif
}

E6t1::~E6t1()
{
    delete ui;
}
