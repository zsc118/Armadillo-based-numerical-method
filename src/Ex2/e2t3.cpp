#include "e2t3.h"
#include "ui_e2t3.h"
#include <QMessageBox>
#include <math.h>
#include <vector>
#include <string>
#include <QFile>
#include <QFileDialog>
#include <QStandardPaths>
#include <QInputDialog>
extern double calStr(std::string s);
namespace _MYFUNCTION {
double distance(const std::vector<double>& x,const std::vector<double>& y)
{
    size_t n(x.size());
    if(!n)return 0;
    if(n!=y.size())return NAN;
    double r(0);
    auto i=x.cend(),j=y.cend();
    do
    {
        double t(*--i-*--j);
        r+=t*=t;
    }while(--n);
    return sqrt(r);
}
}
namespace Ex1 {
extern bool iterative_method(std::vector<double>&,const std::vector<std::string>&, const std::vector<std::string>&);
}
namespace Ex2 {
/*
 * 迭代法
 * current_vec: 迭代向量
 * func       : 迭代函数
 * x          : 未知数
 * epsilon    : 迭代终止条件
 * mid        : 迭代中间结果导出
 *
 * 返回(bool):
 *  true : 迭代失败
 *  false: 迭代成功
 */
bool iterative_method(std::vector<double>& current_vec,const std::vector<std::string> &func, const std::vector<std::string> &x,double epsilon,QString& mid)
{
    for(auto& i:x)mid.append(QString::fromStdString(i)).append(',');
    *(mid.end()-1)='\n';
    for(auto& i:current_vec)mid.append(QString::number(i,'g',14)).append(',');
    std::vector<double> x0(current_vec),x1(x0);
    if(Ex1::iterative_method(current_vec,func,x))
    {
        current_vec=x0;
        return true;
    }
    *(mid.end()-1)='\n';
    for(auto& i:current_vec)mid.append(QString::number(i,'g',14)).append(',');
    while(_MYFUNCTION::distance(current_vec,x1)>=epsilon)
    {
        x1=current_vec;
        if(Ex1::iterative_method(current_vec,func,x))
        {
            current_vec=x0;
            return true;
        }
        *(mid.end()-1)='\n';
        for(auto& i:current_vec)mid.append(QString::number(i,'g',14)).append(',');
    }
    mid.chop(1);
    return false;
}
}
E2t3::E2t3(QWidget *parent) :
    QWidget(parent),ui(new Ui::E2t3),n(0)
{
    ui->setupUi(this);
    ui->tableWidget->setColumnCount(1);
    ui->tableWidget_2->setColumnCount(1);
    ui->tableWidget_3->setColumnCount(1);
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        ui->tableWidget->insertRow(n);
        ui->tableWidget_2->insertRow(n);
        ui->tableWidget_3->insertRow(n++);
    });
    connect(ui->pushButton_2,&QPushButton::clicked,[=](){
        if(!n)QMessageBox::critical(parent,"错误","未知数个数为零！");
        ui->tableWidget->removeRow(--n);
        ui->tableWidget_2->removeRow(n);
        ui->tableWidget_3->removeRow(n);
    });
    connect(ui->pushButton_3,&QPushButton::clicked,[=](){
        if(!n)
        {
            QMessageBox::warning(this,"迭代法","未知数个数为零！");
            return;
        }
        std::vector<std::string> x,f;
        std::vector<double> vec;
        for(unsigned k(0);k<n;++k)
        {
            f.push_back(ui->tableWidget->item(k,0)->text().toStdString());
            x.push_back(ui->tableWidget_2->item(k,0)->text().toStdString());
            double t(calStr(ui->tableWidget_3->item(k,0)->text().toStdString()));
            if(isnan(t))
            {
                QMessageBox::critical(this,"错误","输入有误！");
                return;
            }
            vec.push_back(t);
        }
        if(Ex1::iterative_method(vec,f,x))
        {
            QMessageBox::critical(this,"错误","输入有误！");
            return;
        }
        auto i=vec.end();
        unsigned k(n);
        do ui->tableWidget_3->item(--k,0)->setText(QString::number(*--i));while(k);
    });
    connect(ui->pushButton_4,&QPushButton::clicked,[=](){
        if(!n)
        {
            QMessageBox::warning(this,"迭代法","未知数个数为零！");
            return;
        }
        double eps(calStr(QInputDialog::getText(this,"迭代法","请输入迭代精度：").toStdString()));
        if(isnan(eps))
        {
            QMessageBox::critical(this,"错误","输入有误！");
            return;
        }
        std::vector<std::string> x,f;
        std::vector<double> vec;
        for(unsigned k(0);k<n;++k)
        {
            f.push_back(ui->tableWidget->item(k,0)->text().toStdString());
            x.push_back(ui->tableWidget_2->item(k,0)->text().toStdString());
            double t(calStr(ui->tableWidget_3->item(k,0)->text().toStdString()));
            if(isnan(t))
            {
                QMessageBox::critical(this,"错误","输入有误！");
                return;
            }
            vec.push_back(t);
        }
        QFile file(QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"));
        if(!file.open(QIODevice::WriteOnly))
        {
            QMessageBox::critical(this,"错误","输入有误！");
            return;
        }
        QString s;
        if(Ex2::iterative_method(vec,f,x,eps,s))
        {
            QMessageBox::critical(this,"错误","输入有误！");
            return;
        }
        file.write(s.toStdString().c_str());
        QMessageBox::information(this,"迭代法","导出成功！");
    });
}

E2t3::~E2t3()
{
    delete ui;
}
