#include "e2t2.h"
#include "ui_e2t2.h"
#include <functional>
#include <math.h>
#include <vector>
#include <QFile>
#include <QStandardPaths>
#include <string>
#include <QPushButton>
#include <QMessageBox>
#include <QFileDialog>
#ifdef ZSC_TIME
#include <time.h>
#include <QDebug>
#endif
using namespace std;
extern double calStr(string s);
namespace _MYFUNCTION {
extern string &substring_replace(string &,const string&,const string&,const vector<string>&);
}
namespace CAL {
extern vector<string> func;
}
namespace Ex2 {
/*
 * 二分法
 * f  : 函数
 * a  : 区间左端点
 * b  : 区间右端点
 * eps: 终止条件
 *
 * 返回(double): 近似解
 */
double dichotomy(const function<double(double)>& f,double a,double b,double eps)
{
#ifdef ZSC_TIME
    clock_t start_time=clock();
#endif
    if(a>b)
    {
        double t=a;
        a=b;
        b=t;
    }
//    bool aGe0(f(a)>0),bGe0(f(b)>0);
//    if(!(aGe0^bGe0))return NAN;
    if(f(a)>0)
    {
        if(f(b)>0)return NAN;
        while(b-a>eps)
        {
            double m((a+b)/2);
            if(isnan(f(m)))return NAN;
            if(f(m)>0)
                if(a==m)
                    goto F1;
                else
                    a=m;
            else
                if(b==m)
                    goto F1;
                else
                    b=m;
        }
    }
    else
    {
        if(f(b)<0)return NAN;
        while(b-a>eps)
        {
            double m((a+b)/2);
            if(isnan(f(m)))return NAN;
            if(f(m)>0)
                if(b==m)
                    goto F1;
                else
                    b=m;
            else
                if(a==m)
                    goto F1;
                else
                    a=m;
        }
    }
F1:
#ifdef ZSC_TIME
    clock_t end_time=clock();
    qDebug()<<"二分法用时:"<<end_time-start_time;
#endif
    return (a+b)/2;
}
/*
 * 试位法
 * f  : 函数
 * a  : 区间左端点
 * b  : 区间右端点
 * eps: 终止条件
 *
 * 返回(double): 近似解
 */
double trialPosition(const function<double(double)>& f,double a,double b,double eps)
{
#ifdef ZSC_TIME
    clock_t start_time=clock();
#endif
    if(a>b)
    {
        double t=a;
        a=b;
        b=t;
    }
    if(f(a)>0)
    {
        if(f(b)>0)return NAN;
        while(b-a>eps)
        {
            double m=b-(b-a)*f(b)/(f(b)-f(a));
            if(isnan(f(m)))return NAN;
            if(f(m)>0)
                if(a==m)
                    goto F2;
                else
                    a=m;
            else
                if(b==m)
                    goto F2;
                else
                    b=m;
        }
    }
    else
    {
        if(f(b)<0)return NAN;
        while(b-a>eps)
        {
            double m=b-(b-a)*f(b)/(f(b)-f(a));
            if(isnan(f(m)))return NAN;
            if(f(m)>0)
                if(b==m)
                    goto F2;
                else
                    b=m;
            else
                if(a==m)
                    goto F2;
                else
                    a=m;
        }
    }
F2:
#ifdef ZSC_TIME
    clock_t end_time=clock();
    qDebug()<<"试位法用时:"<<end_time-start_time;
#endif
    return b-(b-a)*f(b)/(f(b)-f(a));
}
/*
 * 二分法
 * f  : 函数
 * a  : 区间左端点
 * b  : 区间右端点
 * eps: 终止条件
 * its: 迭代过程
 *
 * 返回(double): 近似解
 */
double dichotomy(const function<double(double)>& f,double a,double b,double eps,vector<double>& its)
{
#ifdef ZSC_TIME
    clock_t start_time=clock();
#endif
    if(a>b)
    {
        double t=a;
        a=b;
        b=t;
    }
    if(f(a)>0)
    {
        if(f(b)>0)return NAN;
        its.clear();
        while(b-a>eps)
        {
            double m((a+b)/2);
            if(isnan(f(m)))return NAN;
            if(f(m)>0)
                if(a==m)
                    goto F3;
                else
                    a=m;
            else
                if(b==m)
                    goto F3;
                else
                    b=m;
            its.push_back(m);
        }
    }
    else
    {
        if(f(b)<0)return NAN;
        its.clear();
        while(b-a>eps)
        {
            double m((a+b)/2);
            if(isnan(f(m)))return NAN;
            if(f(m)>0)
                if(b==m)
                    goto F3;
                else
                    b=m;
            else
                if(a==m)
                    goto F3;
                else
                    a=m;
            its.push_back(m);
        }
    }
F3:
#ifdef ZSC_TIME
    clock_t end_time=clock();
    qDebug()<<"二分法用时:"<<end_time-start_time;
#endif
    return (a+b)/2;
}
/*
 * 试位法
 * f  : 函数
 * a  : 区间左端点
 * b  : 区间右端点
 * eps: 终止条件
 * its: 迭代过程
 *
 * 返回(double): 近似解
 */
double trialPosition(const function<double(double)>& f,double a,double b,double eps,vector<double>& its)
{
#ifdef ZSC_TIME
    clock_t start_time=clock();
#endif
    if(a>b)
    {
        double t=a;
        a=b;
        b=t;
    }
    if(f(a)>0)
    {
        if(f(b)>0)return NAN;
        its.clear();
        while(b-a>eps)
        {
            double m=b-(b-a)*f(b)/(f(b)-f(a));
            if(isnan(f(m)))return NAN;
            if(f(m)>0)
                if(a==m)
                    goto F4;
                else
                    a=m;
            else
                if(b==m)
                    goto F4;
                else
                    b=m;
            its.push_back(m);
        }
    }
    else
    {
        if(f(b)<0)return NAN;
        its.clear();
        while(b-a>eps)
        {
            double m=b-(b-a)*f(b)/(f(b)-f(a));
            if(isnan(f(m)))return NAN;
            if(f(m)>0)
                if(b==m)
                    goto F4;
                else
                    b=m;
            else
                if(a==m)
                    goto F4;
                else
                    a=m;
            its.push_back(m);
        }
    }
F4:
#ifdef ZSC_TIME
    clock_t end_time=clock();
    qDebug()<<"试位法用时:"<<end_time-start_time;
#endif
    return b-(b-a)*f(b)/(f(b)-f(a));
}
/*
 * 二分法
 * f   : 函数
 * a   : 区间左端点
 * b   : 区间右端点
 * eps : 终止条件
 * its : 迭代过程
 * path: 文件保存路径
 *
 * 返回(double): 近似解
 */
double dichotomy(const function<double(double)>& f,double a,double b,double eps,const QString& path)
{
    QFile file(path);
    if(!file.open(QIODevice::WriteOnly))return NAN;
    if(a>b)
    {
        double t=a;
        a=b;
        b=t;
    }
    file.write(("a,b\n"+to_string(a)+','+to_string(b)).c_str());
    if(f(a)>0)
    {
        if(f(b)>0||isnan(f(b)))return NAN;
        while(b-a>eps)
        {
            double m((a+b)/2),fm(f(m));
            if(isnan(fm))goto F5;
            if(fm>0)
                file.write(('\n'+to_string(a=m)+','+to_string(b)).c_str());
            else
                file.write(('\n'+to_string(a)+','+to_string(b=m)).c_str());
        }
    }
    else
    {
        if(f(b)<0||isnan(f(a))||isnan(f(b)))return NAN;
        while(b-a>eps)
        {
            double m((a+b)/2),fm(f(m));
            if(a==m||b==m||isnan(fm))goto F5;
            if(fm>0)
                file.write(('\n'+to_string(a)+','+to_string(b=m)).c_str());
            else
                file.write(('\n'+to_string(a=m)+','+to_string(b)).c_str());
        }
    }
F5:
    file.close();
    return (a+b)/2;
}
/*
 * 试位法
 * f   : 函数
 * a   : 区间左端点
 * b   : 区间右端点
 * eps : 终止条件
 * its : 迭代过程
 * path: 文件保存路径
 *
 * 返回(double): 近似解
 */
double trialPosition(const function<double(double)>& f,double a,double b,double eps,const QString& path)
{
    QFile file(path);
    if(!file.open(QIODevice::WriteOnly))return NAN;
    if(a>b)
    {
        double t=a;
        a=b;
        b=t;
    }
    file.write(("a,b\n"+to_string(a)+','+to_string(b)).c_str());
    if(f(a)>0)
    {
        if(f(b)>0||isnan(f(b)))return NAN;
        while(b-a>eps)
        {
            double m=b-(b-a)*f(b)/(f(b)-f(a)),fm(f(m));
            if(a==m||b==m||isnan(fm))goto F6;
            if(fm>0)
                file.write(('\n'+to_string(a=m)+','+to_string(b)).c_str());
            else
                file.write(('\n'+to_string(a)+','+to_string(b=m)).c_str());
        }
    }
    else
    {
        if(f(b)<0||isnan(f(a))||isnan(f(b)))return NAN;
        while(b-a>eps)
        {
            double m=b-(b-a)*f(b)/(f(b)-f(a)),fm(f(m));
            if(a==m||b==m||isnan(fm))goto F6;
            if(fm>0)
                file.write(('\n'+to_string(a)+','+to_string(b=m)).c_str());
            else
                file.write(('\n'+to_string(a=m)+','+to_string(b)).c_str());
        }
    }
F6:
    file.close();
    return b-(b-a)*f(b)/(f(b)-f(a));
}
}

E2t2::E2t2(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E2t2)
{
    ui->setupUi(this);
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        string s=ui->plainTextEdit->toPlainText().toStdString();
        double a=calStr(ui->plainTextEdit_2->toPlainText().toStdString()),b=calStr(ui->plainTextEdit_3->toPlainText().toStdString()),eps=calStr(ui->plainTextEdit_4->toPlainText().toStdString());
        if(isnan(a)||isnan(b))
        {
            QMessageBox::critical(this,"错误","输入有误！");
            return;
        }
        auto f=[=](double x)->double{
            string t(s);
            return calStr(_MYFUNCTION::substring_replace(t,"x",to_string(x),CAL::func));
        };
        // 注: 由于to_string函数的限制, 暂时无法计算比1e-6更高的精度
        // 后续可以利用double的内部存储结构(指数位、数字位)将程序改进为自定义转换从而获得更高精度
        if(eps<1e-6)eps=1e-6;
        if(ui->radioButton->isChecked())
            if(isnan(Ex2::dichotomy(f,a,b,eps,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"))))
                QMessageBox::critical(this,"错误","输入有误！");
            else
                QMessageBox::information(this,"图解法","导出成功！");
        else
            if(isnan(Ex2::trialPosition(f,a,b,eps,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"))))
                QMessageBox::critical(this,"错误","输入有误！");
            else
                QMessageBox::information(this,"图解法","导出成功！");
#ifdef ZSC_TIME
        clock_t start_time=clock();
        for(unsigned i(0);i<100000;++i)Ex2::dichotomy(f,a,b,eps,"D:/1.csv");
        clock_t end_time=clock();
        qDebug()<<"二分法耗时:"<<(end_time-start_time)/100000.0;
        start_time=clock();
        for(unsigned i(0);i<100000;++i)Ex2::dichotomy(f,a,b,eps,"D:/1.csv");
        end_time=clock();
        qDebug()<<"试位法耗时:"<<(end_time-start_time)/100000.0;
#endif
    });
}

E2t2::~E2t2()
{
    delete ui;
}
