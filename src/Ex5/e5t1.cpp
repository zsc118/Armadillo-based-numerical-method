#include "e5t1.h"
#include "ui_e5t1.h"
#include <string>
#include <functional>
#include <math.h>
#include <QMessageBox>
using namespace std;
extern double calStr_x(string,const double&);
extern double calStr(string);
#ifdef DEBUG51
#include <QDebug>
#endif
namespace Ex5 {
/*
 * 梯形公式
 * f:被积函数
 * a:积分下限
 * b:积分上限
 */
double trapezoid(function<double(double)> f,const double& a,const double& b)
{
    if(a>=b)throw "区间上限小于区间下限！";
    double fa(f(a));
    if(isnan(fa))throw "积分下限处函数值无法计算！";
    double fb(f(b));
    if(isnan(fb))throw "积分上限处函数值无法计算！";
    return (b-a)*(fa+fb)/2;
}
/*
 * Simpson公式
 * f:被积函数
 * a:积分下限
 * b:积分上限
 */
double Simpson(function<double(double)> f,const double& a,const double& b)
{
    if(a>=b)throw "区间上限小于区间下限！";
    double fa(f(a));
    if(isnan(fa))throw "积分下限处函数值无法计算！";
    double fb(f(b));
    if(isnan(fb))throw "积分上限处函数值无法计算！";
    double fm(f((a+b)/2));
    if(isnan(fm))throw "积分区间中点处函数值无法计算！";
    return (b-a)*(fa+fb+4*fm)/6;
}
/*
 * Newton公式
 * f:被积函数
 * a:积分下限
 * b:积分上限
 */
double Newton(function<double(double)> f,const double& a,const double& b)
{
    if(a>=b)throw "区间上限小于区间下限！";
    double fa(f(a));
    if(isnan(fa))throw "积分下限处函数值无法计算！";
    double fb(f(b));
    if(isnan(fb))throw "积分上限处函数值无法计算！";
    double d0(b-a),d(d0/3),f1(f(a+d));
    if(isnan(f1))throw "积分区间左三等分点处函数值无法计算！";
    double f2(f(b-d));
    if(isnan(f2))throw "积分区间右三等分点处函数值无法计算！";
    return d0*(fa+fb+(f1+f2)*3)/8;
}
/*
 * Cotes公式
 * f:被积函数
 * a:积分下限
 * b:积分上限
 */
double Cotes(function<double(double)> f,const double& a,const double& b)
{
    if(a>=b)throw "区间上限小于区间下限！";
    double fa(f(a));
    if(isnan(fa))throw "积分下限处函数值无法计算！";
    double fb(f(b));
    if(isnan(fb))throw "积分上限处函数值无法计算！";
    double d0(b-a),d(d0/4),x(a+d),f1(f(x));
    if(isnan(f1))throw "积分区间左四等分点处函数值无法计算！";
    double fm(f(x+=d));
    if(isnan(fm))throw "积分区间中点处函数值无法计算！";
    double f2(f(x+=d));
    if(isnan(f2))throw "积分区间右四等分点处函数值无法计算！";
    return d0*(7*(fa+fb)+32*(f1+f2)+12*fm)/90;
}
}
E5t1::E5t1(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E5t1)
{
    ui->setupUi(this);
#ifdef DEBUG51
    double pi(52163.0 / 16604),pi4(pi/4),ans(1-sqrt(2)/2);
    double (*f)(double)=sin;
    try {
        double x1(Ex5::trapezoid(f,0.0,pi4)),e1(x1-ans),x2(Ex5::Simpson(f,0.,pi4)),e2(x2-ans),x3(Ex5::Newton(f,0.,pi4)),e3(x3-ans),x4(Ex5::Cotes(f,0.,pi4)),e4(x4-ans);
        qDebug()<<x1<<'\t'<<e1<<'\n'<<x2<<'\t'<<e2<<'\n'<<x3<<'\t'<<e3<<'\n'<<x4<<'\t'<<e4<<'\n'<<ans;
    } catch (const char* s) {
        qDebug()<<s;
    }
#else
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        try {
            double a(calStr(ui->plainTextEdit_2->toPlainText().toStdString()));
            if(isnan(a))throw "积分下限输入有误！";
            double b(calStr(ui->plainTextEdit_3->toPlainText().toStdString()));
            if(isnan(b))throw "积分上限输入有误！";
            string s(ui->plainTextEdit->toPlainText().toStdString());
            auto f=[=](double x)->double{
                return calStr_x(s,x);
            };
            if(ui->radioButton->isChecked())
                ui->label_4->setText("I="+QString::number(Ex5::trapezoid(f,a,b)));
            else if(ui->radioButton_2->isChecked())
                ui->label_4->setText("I="+QString::number(Ex5::Simpson(f,a,b)));
            else if(ui->radioButton_3->isChecked())
                ui->label_4->setText("I="+QString::number(Ex5::Newton(f,a,b)));
            else
                ui->label_4->setText("I="+QString::number(Ex5::Cotes(f,a,b)));
        } catch (const char* s) {
            QMessageBox::critical(this,"插值型求积",s);
        }
    });
#endif
}

E5t1::~E5t1()
{
    delete ui;
}
