#include "e2t1.h"
#include "ui_e2t1.h"
#include <QStandardPaths>
#include <functional>
#include <QFile>
#include <math.h>
#include <QMessageBox>
#include <QFileDialog>
#include <string>
#ifdef ZSC_TIME
#include <time.h>
#include <QDebug>
#endif
using namespace std;
namespace CAL {
extern std::vector<std::string> func;
}
extern double calStr(string s);
namespace _MYFUNCTION {
extern string &substring_replace(string &,const string&,const string&,const vector<string>&);
}
namespace Ex2 {
/*
 * 导出画图数据(一元函数)
 * f     : 函数
 * a     : 画图区间左端点
 * b     : 画图区间右端点
 * n     : 区间划分数
 * path  : 文件保存路径
 * title : 表头
 *
 * 返回(bool):
 *  true  :导出失败
 *  false :导出成功
 */
bool plot(function<double(double)> f,double a=0,double b=1,unsigned n=1000,QString path=QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),const char*title="x,f(x)")
{
#ifdef ZSC_TIME
    clock_t start_time=clock();
#endif
    QFile file(path);
    if (!file.open(QIODevice::WriteOnly))return true;
    if(b<a)
    {
        double t(a);
        a=b;
        b=t;
    }
    file.write(title);
    double step((b-a)/n++);
    do
    {
        double f_a(f(a));
        if(isnan(f_a))
        {
            file.close();
            return true;
        }
        file.write(('\n'+QString::number(a)+','+QString::number(f_a)).toUtf8());
        a+=step;
    }while(--n);
#ifdef ZSC_TIME
    clock_t end_time=clock();
    qDebug()<<"导出画图数据(一元函数)用时:"<<end_time-start_time;
#endif
    return false;
}
/*
 * 导出画图数据(一元函数)
 * f     : 函数
 * a     : 画图区间左端点
 * b     : 画图区间右端点
 * step  : 步长
 * path  : 文件保存路径
 * title : 表头
 *
 * 返回(bool):
 *  true  :导出失败
 *  false :导出成功
 */
bool plot(function<double(double)> f,double a,double b,double step,QString path=QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),const char*title="x,f(x)")
{
#ifdef ZSC_TIME
    clock_t start_time=clock();
#endif
    QFile file(path);
    if (!file.open(QIODevice::WriteOnly))return true;
    if(b<a)
    {
        double t(a);
        a=b;
        b=t;
    }
    file.write(title);
    do
    {
        double f_a(f(a));
        if(isnan(f_a))
        {
            file.close();
            return true;
        }
        file.write(('\n'+QString::number(a)+','+QString::number(f_a)).toUtf8());
    }while((a+=step)<=b);
#ifdef ZSC_TIME
    clock_t end_time=clock();
    qDebug()<<"导出画图数据(一元函数)用时:"<<end_time-start_time;
#endif
    return false;
}
}
E2t1::E2t1(QWidget *parent) :
    QWidget(parent), ui(new Ui::E2t1)
{
    ui->setupUi(this);
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        double a(calStr(ui->plainTextEdit_2->toPlainText().toStdString())),b(calStr(ui->plainTextEdit_3->toPlainText().toStdString()));
        unsigned n(ui->plainTextEdit_4->toPlainText().toUInt());
        if(!n||isnan(a)||isnan(b))
        {
            QMessageBox::critical(this,"错误","输入有误！");
            return;
        }
        string s_f(ui->plainTextEdit->toPlainText().toStdString());
        auto f=[=](double x)->double{
            string s(s_f);
            return calStr(_MYFUNCTION::substring_replace(s,"x",to_string(x),CAL::func));
        };
        if(Ex2::plot(f,a,b,n,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()")))QMessageBox::critical(this,"错误","输入有误！");else QMessageBox::information(this,"图解法","导出成功！");
    });
}

E2t1::~E2t1()
{
    delete ui;
}
