#include "e5t2.h"
#include "ui_e5t2.h"
#include <functional>
#include <string>
#include <math.h>
#include <QMessageBox>
#include <QInputDialog>
using namespace std;
extern double calStr_x(string,const double&);
extern double calStr(string);
#ifdef DEBUG52
#include <QDebug>
#endif
#ifdef TIME52
#include <QDebug>
#include <time.h>
#endif
namespace Ex5 {
/*
 * 复化梯形公式
 * f :被积函数
 * a :积分下限
 * b :积分上限
 * e :精度
 */
double composite_trapezoid(function<double(double)> f,const double& a,const double& b,const double& e=1e-6)
{
    if(a>=b)throw "区间上限小于区间下限！";
    double I(f(a));
    if(isnan(I)||isnan(I+=f(b)))throw "区间内存在奇点！";
    double step(b-a),pre_I(I*=step/=2);
    if(isnan((I/=2)+=step*f(a+step)))throw "区间内存在奇点！";
    while(abs(I-pre_I)>e)
    {
        double h(step),x(a+(step/=2)),s(f(x));
        if(isnan(s))throw "区间内存在奇点！";
        while((x+=h)<b)
            if(isnan(s+=f(x)))
                throw "区间内存在奇点！";
        pre_I=I;
        (I/=2)+=step*s;
    }
    return I;
}
/*
 * 复化梯形公式(端点可能为奇点)
 * f :被积函数
 * a :积分下限
 * b :积分上限
 * fa:积分下限处的函数值
 * fb:积分上限处的函数值
 * e :精度
 */
double composite_trapezoid(function<double(double)> f,const double& a,const double& b,const double& fa,const double& fb,const double& e=1e-6)
{
    if(a>=b)throw "区间上限小于区间下限！";
    double step(b-a),I(fa+fb),pre_I(I*=step/=2);
    if(isnan((I/=2)+=step*f(a+step)))throw "区间内存在奇点！";
    while(abs(I-pre_I)>e)
    {
        double h(step),x(a+(step/=2)),s(f(x));
        if(isnan(s))throw "区间内存在奇点！";
        while((x+=h)<b)
            if(isnan(s+=f(x)))
                throw "区间内存在奇点！";
        pre_I=I;
        (I/=2)+=step*s;
    }
    return I;
}
/*
 * 复化Simpson公式(端点可能为奇点)
 * f :被积函数
 * a :积分下限
 * b :积分上限
 * fa:积分下限处的函数值
 * fb:积分上限处的函数值
 * e :精度
 */
double composite_Simpson(function<double(double)> f,const double& a,const double& b,const double& fa,const double& fb,const double& e=1e-6)
{
    if(a>=b)throw "区间上限小于区间下限！";
    double step(b-a),m(a+(step/=2)),fx(f(m));
    if(isnan(fx))throw "区间内存在奇点！";
    double I((fx*4+fa+fb)/6),pre_I(I),pre_fx(fx),h(step);
    fx=f(a+(step/=2));
    if(isnan((I/=2)+=((fx+=f(b-step))*2-pre_fx)*h/3))throw "区间内存在奇点！";
    while(abs(I-pre_I)>e)
    {
        pre_fx=fx,pre_I=I,h=step;
        if(isnan(fx=f(m=a+(step/=2))))throw "区间内存在奇点！";
        while((m+=h)<b)
            if(isnan(fx+=f(m)))
                throw "区间内存在奇点！";
        (I/=2)+=(fx*2-pre_fx)*h/3;
    }
    return I;
}
/*
 * 复化Simpson公式
 * f :被积函数
 * a :积分下限
 * b :积分上限
 * e :精度
 */
inline double composite_Simpson(function<double(double)> f,const double& a,const double& b,const double& e=1e-6)
{
    double fa(f(a));
    if(isnan(fa))throw "区间内存在奇点！";
    double fb(f(b));
    if(isnan(fb))throw "区间内存在奇点！";
    return composite_Simpson(f,a,b,fa,fb,e);
}
}
E5t2::E5t2(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E5t2)
{
    ui->setupUi(this);
#ifdef DEBUG52
    auto f=[](double x)->double{return sin(x)/x;};
    try {
        double I1(Ex5::composite_trapezoid(f,0.,1.,1.,sin(1.),1e-8));
        qDebug()<<I1;
        double I2(Ex5::composite_Simpson(f,0.,1.,1.,sin(1.),1e-8));
        qDebug()<<I2;
    } catch (const char* s) {
        qDebug()<<s;
    }
#else
#ifdef TIME52
    auto f=[](double x)->double{return sin(x)/x;};
    const unsigned short end(-1);
    {
        auto start_time=clock();
        for(unsigned short i(0);i!=end;++i)
            Ex5::composite_trapezoid(f,0.,1.,1.,sin(1.),1e-8);
        auto end_time=clock();
        qDebug()<<"复化梯形耗时:"<<(end_time-start_time)/65536.0;
    }
    {
        auto start_time=clock();
        for(unsigned short i(0);i!=end;++i)
            Ex5::composite_Simpson(f,0.,1.,1.,sin(1.),1e-8);
        auto end_time=clock();
        qDebug()<<"复化Simpson耗时:"<<(end_time-start_time)/65536.0;
    }
#endif
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        try {
            double a(calStr(ui->plainTextEdit_2->toPlainText().toStdString()));
            if(isnan(a))throw "积分下限输入有误！";
            double b(calStr(ui->plainTextEdit_3->toPlainText().toStdString()));
            if(isnan(b))throw "积分下限输入有误！";
            double e(calStr(ui->plainTextEdit_4->toPlainText().toStdString()));
            if(isnan(e))throw "求积精度输入有误！";
            string s(ui->plainTextEdit->toPlainText().toStdString());
            auto f=[s](double x)->double{
                return calStr_x(s,x);
            };
            double fa(f(a));
            if(isnan(fa))
                if(isnan(fa=calStr(QInputDialog::getText(this,"复化求积法","积分下限处为奇点，请输入积分下限处的函数值:").toStdString())))
                    throw "函数值输入有误！";
            double fb(f(b));
            if(isnan(fb))
                if(isnan(fb=calStr(QInputDialog::getText(this,"复化求积法","积分上限处为奇点，请输入积分上限处的函数值:").toStdString())))
                    throw "函数值输入有误！";
            if(ui->radioButton->isChecked())
                QMessageBox::information(this,"复化求积法","积分近似值为:"+QString::number(Ex5::composite_trapezoid(f,a,b,fa,fb,e)));
            else
                QMessageBox::information(this,"复化求积法","积分近似值为:"+QString::number(Ex5::composite_Simpson(f,a,b,fa,fb,e)));
        } catch (const char* s) {
            QMessageBox::critical(this,"复化求积法",s);
        }
    });
#endif
}

E5t2::~E5t2()
{
    delete ui;
}
