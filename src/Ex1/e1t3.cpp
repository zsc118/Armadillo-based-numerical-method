#include "e1t3.h"
#include "ui_e1t3.h"
#include <QString>
#include <armadillo>
#include <QStandardPaths>
#include <QFile>
#include <string>
#include <QFileDialog>
#include <QMessageBox>
using namespace arma;
using namespace std;
namespace Ex1 {
/*
 * 线性随机IFS迭代
 * ifs : IFS码
 * n   : 迭代次数
 * x0  : 迭代初值
 * path: 迭代过程保存路径
 * e   : 实数比较精度
 *
 * 返回(bool):
 *  true : 迭代失败
 *  false: 迭代成功
 */
bool IFS_iteration(const mat& ifs,unsigned n,const QString& path=QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),vec x0=randu(2),const double &e=1e-6)
{
    if(ifs.n_cols!=7||x0.n_elem!=2)return true;
    if(!n)return true;
    vec p=cumsum(ifs.col(6));
    if(abs(p.at(p.n_elem-1)-1.)>e)return true;
    QFile file(path);
    if(!file.open(QIODevice::WriteOnly))return true;
    file.write(("x,y\n"+to_string(x0.at(0))+','+to_string(x0.at(1))).c_str());
    static default_random_engine engine(time(nullptr));
    static uniform_real_distribution<double> distribution;
    do
    {
        double d(distribution(engine));
        unsigned i(0);
        while(p.at(i)<=d)
            if(++i==p.n_elem)
            {
                --i;
                break;
            }
        d=ifs.at(i,0)*x0.at(0)+ifs.at(i,1)*x0.at(1)+ifs.at(i,4);
        x0.at(1)=ifs.at(i,2)*x0.at(0)+ifs.at(i,3)*x0.at(1)+ifs.at(i,5);
        x0.at(0)=d;
        file.write(('\n'+to_string(x0.at(0))+','+to_string(x0.at(1))).c_str());
    }while(--n);
    file.close();
    return false;
}
}

E1t3::E1t3(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E1t3)
{
    ui->setupUi(this);
    ui->widget->setTitle("IFS码");
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        if(Ex1::IFS_iteration(ui->widget->x,ui->plainTextEdit->toPlainText().toUInt(),QFileDialog::getSaveFileName(this,"迭代过程保存",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()")))QMessageBox::critical(this,"错误","输入错误！");else QMessageBox::information(this,"线性随机IFS迭代法","导出成功！");
    });
}

E1t3::~E1t3()
{
    delete ui;
}
