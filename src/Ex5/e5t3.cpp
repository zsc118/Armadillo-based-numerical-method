#include "e5t3.h"
#include "ui_e5t3.h"
#include <functional>
#include <math.h>
#include <string>
#include <QMessageBox>
#include <QInputDialog>
using namespace std;
typedef function<double(double)> func;
extern double calStr_x(string,const double&);
extern double calStr(string);
#ifdef DEBUG53
#include <QDebug>
#endif
namespace Ex5 {
class Romberg_Table
{
    double *arr;
    unsigned n,k;
public:
    Romberg_Table(const func& f,const double& a,const double& b,const double& fa,const double& fb,unsigned n0):n(n0),k(0)
    {
        if(!n0)throw "n不能为零！";
        double* p(arr=new double[n0]),h(b-a);
        *p=h*(fa+fb)/2;
        unsigned i(n0);
        while(--i)
        {
            double pre_h(h),x(a+(h/=2)),fx(f(x)),&I(*p);
            if(isnan(fx))
            {
                delete[]arr;
                throw "区间内存在奇点！";
            }
            while((x+=pre_h)<b)
                if(isnan(fx+=f(x)))
                {
                    delete[]arr;
                    throw "区间内存在奇点！";
                }
            *++p=I/2+fx*h;
        }
    }
    ~Romberg_Table()
    {
        delete[] arr;
    }
    //外推, 注意n不能等于1
    void extrapolation()noexcept
    {
        unsigned k4(1<<(++k<<1)),k41(k4-1),m(--n);//k4表示4^k
        double*p(arr);
        *p=(k4**(p+1)-*p)/k41;
        while(--m)
        {
            ++p;
            *p=(k4**(p+1)-*p)/k41;
        }
    }
    //输出当前值
    inline double& get_value()noexcept
    {
        return arr[n-1];
    }
    //判断是否达到精度要求
    inline bool judge_accuracy(const double& e)noexcept
    {
        if(n==1)return false;
        return abs(arr[n-2]-arr[n-1])>=e;
    }
#ifdef DEBUG53
    //输出表格
    void print_table()
    {
        for(unsigned i(0);i!=n;++i)
            qDebug()<<arr[i];
        qDebug();
    }
#endif
};
/*
 * Romberg求积法(端点可能为奇点)
 * f :被积函数
 * a :积分下限
 * b :积分上限
 * fa:积分下限处的函数值
 * fb:积分上限处的函数值
 * e :精度
 * n :复化次数
 */
double Romberg(func f,const double& a,const double& b,const double& fa,const double& fb,const double& e=1e-8,const unsigned& n=16)
{
    if(a>=b)throw "区间上限小于区间下限！";
    if(!n)throw "复化次数输入有误！";
    Romberg_Table table(f,a,b,fa,fb,n);
#ifdef DEBUG53
    table.print_table();
#endif
    while(table.judge_accuracy(e))
    {
        table.extrapolation();
#ifdef DEBUG53
        table.print_table();
#endif
    }
    return table.get_value();
}
/*
 * Romberg求积法
 * f :被积函数
 * a :积分下限
 * b :积分上限
 * e :精度
 * n :复化次数
 */
double Romberg(func f,const double& a,const double& b,const double& e=1e-8,const unsigned& n=16)
{
    double fa(f(a));
    if(isnan(fa))throw "积分下限为奇点！";
    double fb(f(b));
    if(isnan(fb))throw "积分上限为奇点！";
    return Romberg(f,a,b,fa,fb,e,n);
}
}
E5t3::E5t3(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E5t3)
{
    ui->setupUi(this);
#ifdef DEBUG53
    auto f=[](double x)->double{
        return sqrt(x)*log(x);
    };
    try {
        qDebug()<<Ex5::Romberg(f,0.,1.,0.,0.,1e-6,15);
    } catch (const char* s) {
        qDebug()<<s;
    }
#else
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        try {
            double a(calStr(ui->lineEdit_2->text().toStdString()));
            if(isnan(a))throw "积分下限输入有误！";
            double b(calStr(ui->lineEdit_3->text().toStdString()));
            if(isnan(b))throw "积分下限输入有误！";
            double e(calStr(ui->lineEdit_4->text().toStdString()));
            if(isnan(e))throw "求积精度输入有误！";
            unsigned n(ui->lineEdit_5->text().toUInt());
            string s(ui->lineEdit->text().toStdString());
            auto f=[s](double x)->double{
                return calStr_x(s,x);
            };
            double fa(f(a));
            if(isnan(fa))
                if(isnan(fa=calStr(QInputDialog::getText(this,"Romberg求积法","积分下限处为奇点，请输入积分下限处的函数值:").toStdString())))
                    throw "函数值输入有误！";
            double fb(f(b));
            if(isnan(fb))
                if(isnan(fb=calStr(QInputDialog::getText(this,"Romberg求积法","积分上限处为奇点，请输入积分上限处的函数值:").toStdString())))
                    throw "函数值输入有误！";
            QMessageBox::information(this,"Romberg求积法","积分近似值:"+QString::number(Ex5::Romberg(f,a,b,fa,fb,e,n)));
        } catch (const char* s) {
            QMessageBox::critical(this,"Romberg求积法",s);
        }
    });
#endif
}

E5t3::~E5t3()
{
    delete ui;
}
