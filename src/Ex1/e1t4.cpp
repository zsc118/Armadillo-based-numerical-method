#include "e1t4.h"
#include "ui_e1t4.h"
#include <armadillo>
#include <QFile>
#include <string>
#include <QMessageBox>
#include <QFileDialog>
#include <QStandardPaths>
using namespace std;
using namespace arma;
extern double calStr(std::string s);
namespace Ex1 {
/*
 * Julia集绘制
 * [xa,xb]: 绘图区域的x范围
 * [ya,yb]: 绘图区域的y范围
 * c      : 初始点
 * M      : 逃逸半径
 * N      : 最大迭代数
 * n      : 区间分割数
 */
umat Julia(double xa,double xb,double ya,double yb,cx_double c=cx_double(-0.46,0.57),double M=2,unsigned N=100,unsigned n=512)
{
    ++n;
    umat R=zeros<umat>(n,n);
    rowvec x(linspace(xa,xb,n));
    vec y(linspace(ya,yb,n));
    cx_mat z(repmat(x,n,1),repmat(y,1,n));
    unsigned k(0);
    uvec ind=find(abs(z)>M);
    for(auto &i:ind)z.at(i)=NAN;
    while(++k!=N)
    {
        ind=find(abs((z%=z)+=c)>M);
        for(auto &i:ind)
        {
            z.at(i)=NAN;
            R.at(i)=k;
        }
    }
    ind=find_finite(z);
    for(auto &i:ind)R.at(i)=N;
    return R;
}
/*
 * Julia集绘制
 * [xa,xb]: 绘图区域的x范围
 * [ya,yb]: 绘图区域的y范围
 * c      : 初始点
 * M      : 逃逸半径
 * N      : 最大迭代数
 * n      : 区间分割数
 * path   : 文件保存路径
 *
 * 返回(bool):
 *  true : 导出失败
 *  false: 导出成功
 */
bool Julia(double xa,double xb,double ya,double yb,cx_double c,double M,unsigned N,unsigned n,const QString& path)
{
    QFile file(path);
    if(!file.open(QIODevice::WriteOnly))return true;
    file.write("x,y,W");
    ++n;
    umat R=zeros<umat>(n,n);
    rowvec x(linspace<rowvec>(xa,xb,n));
    vec y(linspace(ya,yb,n));
    cx_mat z(repmat(x,n,1),repmat(y,1,n));
    unsigned k(0);
    uvec ind=find(abs(z)>M);
    for(auto &i:ind)z.at(i)=NAN;
    while(++k!=N)
    {
        ind=find(abs((z%=z)+=c)>M);
        for(auto &i:ind)
        {
            z.at(i)=NAN;
            R.at(i)=k;
        }
    }
    ind=find_finite(z);
    for(auto &i:ind)R.at(i)=N;
    k=0;
    do
    {
        N=0;
        string s('\n'+to_string(x.at(k))+',');
        do file.write((s+to_string(y.at(N))+','+to_string(R.at(k,N))).c_str());while(++N!=n);
    }while(++k!=n);
    file.close();
    return false;
}
/*
 * Mandelbrot集绘制
 * [xa,xb]: 绘图区域的x范围
 * [ya,yb]: 绘图区域的y范围
 * M      : 逃逸半径
 * N      : 最大迭代数
 * n      : 区间分割数
 */
umat Mandelbrot(double xa,double xb,double ya,double yb,double M=2,unsigned N=100,unsigned n=512)
{
    ++n;
    umat R=zeros<umat>(n,n);
    rowvec x(linspace(xa,xb,n));
    vec y(linspace(ya,yb,n));
    cx_mat z(repmat(x,n,1),repmat(y,1,n)),c(z);
    unsigned k(0);
    uvec ind=find(abs(z)>M);
    for(auto &i:ind)z.at(i)=NAN;
    while(++k!=N)
    {
        ind=find(abs((z%=z)+=c)>M);
        for(auto &i:ind)
        {
            z.at(i)=NAN;
            R.at(i)=k;
        }
    }
    ind=find_finite(z);
    for(auto &i:ind)R.at(i)=N;
    return R;
}

bool Mandelbrot(double xa,double xb,double ya,double yb,double M,unsigned N,unsigned n,const QString& path)
{
    QFile file(path);
    if(!file.open(QIODevice::WriteOnly))return true;
    file.write("x,y,W");
    ++n;
    umat R=zeros<umat>(n,n);
    rowvec x(linspace<rowvec>(xa,xb,n));
    vec y(linspace(ya,yb,n));
    cx_mat z(repmat(x,n,1),repmat(y,1,n)),c(z);
    unsigned k(0);
    uvec ind=find(abs(z)>M);
    for(auto &i:ind)z.at(i)=NAN;
    while(++k!=N)
    {
        ind=find(abs((z%=z)+=c)>M);
        for(auto &i:ind)
        {
            z.at(i)=NAN;
            R.at(i)=k;
        }
    }
    ind=find_finite(z);
    for(auto &i:ind)R.at(i)=N;
    k=0;
    do
    {
        N=0;
        string s('\n'+to_string(x.at(k))+',');
        do file.write((s+to_string(y.at(N))+','+to_string(R.at(k,N))).c_str());while(++N!=n);
    }while(++k!=n);
    file.close();
    return false;
}
}

E1t4::E1t4(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E1t4)
{
    ui->setupUi(this);
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        double xa=calStr(ui->plainTextEdit->toPlainText().toStdString()),xb=calStr(ui->plainTextEdit_2->toPlainText().toStdString()),ya=calStr(ui->plainTextEdit_4->toPlainText().toStdString()),yb=calStr(ui->plainTextEdit_3->toPlainText().toStdString()),M=calStr(ui->plainTextEdit_7->toPlainText().toStdString());
        cx_double c(calStr(ui->plainTextEdit_6->toPlainText().toStdString()),calStr(ui->plainTextEdit_5->toPlainText().toStdString()));
        unsigned N=ui->plainTextEdit_8->toPlainText().toUInt(),n=ui->plainTextEdit_9->toPlainText().toUInt();
        if(!N||!n||isnan(xa)||isnan(xb)||isnan(ya)||isnan(yb)||isnan(c.real()||isnan(c.imag()))||isnan(M))
        {
            QMessageBox::critical(this,"错误","输入有误！");
            return;
        }
        if(ui->radioButton->isChecked())
            if(Ex1::Julia(xa,xb,ya,yb,c,M,N,n,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()")))
                QMessageBox::critical(this,"错误","输入有误！");
            else
                QMessageBox::information(this,"复平面上的迭代","导出成功！");
        else
            if(Ex1::Mandelbrot(xa,xb,ya,yb,M,N,n,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()")))
                QMessageBox::critical(this,"错误","输入有误！");
            else
                QMessageBox::information(this,"复平面上的迭代","导出成功！");
    });
}

E1t4::~E1t4()
{
    delete ui;
}
