#include "e4t1.h"
#include "ui_e4t1.h"
#include <vector>
#include <unordered_set>
#include <QFile>
#include <QMessageBox>
#include <QFileDialog>
#include <QStandardPaths>
#include <string>
using namespace arma;
using namespace std;
#ifdef DEBUG41
#include <QDebug>
#endif
#ifdef TIME41
#include <QDebug>
#include <time.h>
#endif
namespace _MYFUNCTION {
/*
 * 检查向量是否包含重复元素
 *  包含  :true
 *  不包含:false
 */
bool check_unique(const vec& x)
{
    unordered_set<double> s;
    for(const auto& i:x)
        if(!s.insert(i).second)
            return true;
    return false;
}
}
namespace Ex4 {
/*
 * Newton插值
 * x0:观测点
 * y0:观测数据
 * x1:插值点
 * y1:预测值
 */
void Newton(const vec& x0,const vec& y0,const vec& x1,vec& y1)
{
    if(x1.empty())
    {
        y1.clear();
        return;
    }
    if(x0.empty()||y0.empty())throw "观测集为空！";
    if(x0.n_elem!=y0.n_elem)throw "观测点与观测数据向量大小不一致！";
    if(x0.n_elem<2)throw "观测点数量小于2！";
    if(_MYFUNCTION::check_unique(x0))throw "观测集包含重复元素！";
    unsigned k(x0.n_elem),l(k),g(0);
    double*a(new double[k]),*p(a+k);
    auto q(y0.cend());
    do *--p=*--q;while(--l);
    y1.set_size(x1.n_elem).fill(*p);
    vec xt(ones(x1.n_elem));
    while(l=--k)
    {
        xt%=x1-x0.at(g);
        p=a;
        q=x0.cbegin();
        y1+=((*p-=*(p+1))/=*q-*(q+(++g)))*xt;
        while(--l)
        {
            ++q,++p;
            (*p-=*(p+1))/=*q-*(q+g);
        }
    }
    delete[]a;
}
/*
 * Lagrange插值
 * x0:观测点
 * y0:观测数据
 * x1:插值点
 * y1:预测值
 */
void Lagrange(const vec& x0,const vec& y0,const vec& x1,vec& y1)
{
    if(x1.empty())
    {
        y1.clear();
        return;
    }
    if(x0.empty()||y0.empty())throw "观测集为空！";
    if(x0.n_elem!=y0.n_elem)throw "观测点与观测数据向量大小不一致！";
    if(x0.n_elem<2)throw "观测点数量小于2！";
    if(_MYFUNCTION::check_unique(x0))throw "观测集包含重复元素！";
    y1.zeros(x1.n_elem);
    const unsigned start(-1);
    for(unsigned i(0);i!=x1.n_elem;++i)
        for(unsigned j(0);j!=x0.n_elem;++j)
        {
            double t(y0.at(j));
            unsigned k(start);
            while(++k!=j)(t*=x1.at(i)-x0.at(k))/=x0.at(j)-x0.at(k);
            while(++k!=x0.n_elem)(t*=x1.at(i)-x0.at(k))/=x0.at(j)-x0.at(k);
            y1.at(i)+=t;
        }
}
}
E4t1::E4t1(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E4t1)
{
    ui->setupUi(this);
#ifdef DEBUG41
    double (*p)(double)=sqrt;
    vec x0{2.0,2.1,2.2},y0(x0),x1(1),y1;
    y0.transform(p);
    x1.at(0)=2.15;
    Ex4::Newton(x0,y0,x1,y1);
    qDebug()<<y1.at(0);
    Ex4::Lagrange(x0,y0,x1,y1);
    qDebug()<<y1.at(0);
#else
#ifdef TIME41
    double (*p)(double)=sqrt;
    vec x0(linspace(2,2.3)),y0(x0),x1(1),y1;
    y0.transform(p);
    x1.at(0)=2.15;
    const unsigned short end(-1);
    {
        auto start_time=clock();
        for(unsigned short i(0);i!=end;++i)
            Ex4::Newton(x0,y0,x1,y1);
        auto end_time=clock();
        qDebug()<<"Newton耗时:"<<(end_time-start_time)/65536.0;
    }
    {
        auto start_time=clock();
        for(unsigned short i(0);i!=end;++i)
            Ex4::Lagrange(x0,y0,x1,y1);
        auto end_time=clock();
        qDebug()<<"Lagrange耗时:"<<(end_time-start_time)/65536.0;
    }
#else
    ui->widget->setTitle("观测点");
    ui->widget_2->setTitle("观测数据");
    ui->widget_3->setTitle("插值点");
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        if(ui->widget_3->x.empty())
        {
            QMessageBox::warning(this,"多项式插值","插值集为空！");
            return;
        }
        vec t;
        try {
            if(ui->radioButton->isChecked())
                Ex4::Lagrange(ui->widget->x,ui->widget_2->x,ui->widget_3->x,t);
            else
                Ex4::Newton(ui->widget->x,ui->widget_2->x,ui->widget_3->x,t);
        } catch (const char* s) {
            QMessageBox::critical(this,"多项式插值",s);
            return;
        }
        QFile file(QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"));
        if(!file.open(QIODevice::WriteOnly))
        {
            QMessageBox::critical(this,"多项式插值","文件保存失败！");
            return;
        }
        file.write("x,y");
        for(unsigned i(0);i!=t.n_elem;++i)
        {
            string s("\n");
            file.write((((s+=to_string(ui->widget_3->x.at(i)))+=',')+=to_string(t.at(i))).c_str());
        }
        file.close();
        QMessageBox::information(this,"多项式插值","导出成功！");
    });
#endif
#endif
}

E4t1::~E4t1()
{
    delete ui;
}
