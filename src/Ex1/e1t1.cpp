#include "e1t1.h"
#include "ui_e1t1.h"
#include <QString>
#include <string>
#include <math.h>
#include <QMessageBox>
extern double calStr(std::string s);
namespace Ex1 {
/*
 * 十进制展开
 * a        : 数字部分
 * e        : 指数部分
 * x        : 要展开的十进制数
 * precision: 精度
 *
 * 注：不考虑x的正负
 */
void decimal_expansion(double x,unsigned char* &a,short&e,int precision=14)
{
    if(x<0)x=-x;
    QString s(QString::number(x,'E',precision));
    unsigned char*p(a=new unsigned char[++precision]);
    e=s.mid(precision+2).toShort()+1;
    QString::iterator q(s.begin());
    *p=(q++)->toLatin1();
    while(--precision)*++p=(++q)->toLatin1();
}
// 只计算指数部分
short decimal_expansion(double x,int precision=14)
{
    if(x<0)x=-x;
    QString s(QString::number(x,'E',precision));
    return s.mid(precision+2).toShort()+1;
}
/*
 * x0     : 真值
 * x      : 测量值
 * N      : 有效数字
 * epsilon: 误差限
 */
void cal_N_epsilon(const double& x0,const double& x,short &N,double &epsilon,int precision=14)
{
    double err(x0-x);
    short m(decimal_expansion(x0,precision)),q;
    unsigned char*a;
    decimal_expansion(err,a,q,precision);
    if(*a<'5')
        N=m-q;
    else
        N=m-q-1;
    epsilon=pow(10,m-N)/2;
    delete a;
}
}

E1t1::E1t1(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E1t1)
{
    ui->setupUi(this);
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        double x0(calStr(ui->plainTextEdit->toPlainText().toStdString())),x(calStr(ui->plainTextEdit_2->toPlainText().toStdString())),epsilon;
        if(isnan(x0)||isnan(x))
        {
            QMessageBox::critical(this,"错误","输入有误！\n请重新输入！");
            return;
        }
        short N;
        Ex1::cal_N_epsilon(x0,x,N,epsilon);
        ui->label_3->setText("有效数字:"+QString::number(N)+"\t\t误差限："+QString::number(epsilon));
    });
}

E1t1::~E1t1()
{
    delete ui;
}
