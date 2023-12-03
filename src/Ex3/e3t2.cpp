#include "e3t2.h"
#include "ui_e3t2.h"
#include <QMessageBox>
using namespace arma;
namespace Ex3 {
/*
 * LU分解
 * L:下三角矩阵
 * U:上三角矩阵
 * A:待分解矩阵
 * e:精度
 *
 * 返回(bool):
 *  true : 分解失败
 *  false: 分解成功
 */
bool LU(mat& L,mat& U,const mat& A,const double& e=1e-6)
{
    if(A.n_cols==1)
    {
        L.ones(1,1);
        U.resize(1,1);
        if(abs(U.at(0,0)=A.at(0,0))<e)return true;
        return false;
    }
    L.eye(A.n_cols,A.n_cols);
    U.zeros(A.n_cols,A.n_cols);
    unsigned n(A.n_cols-1);
    for(unsigned i(0);i!=n;++i)
    {
        U.at(i,i)=A.at(i,i);
        for(unsigned k(0);k!=i;++k)
            U.at(i,i)-=L.at(i,k)*U.at(k,i);
        if(abs(U.at(i,i))<e)return true;
        for(unsigned j(i+1);j!=A.n_cols;++j)
        {
            L.at(j,i)=A.at(j,i);
            U.at(i,j)=A.at(i,j);
            for(unsigned k(0);k!=i;++k)
            {
                U.at(i,j)-=L.at(i,k)*U.at(k,j);
                L.at(j,i)-=L.at(j,k)*U.at(k,i);
            }
            L.at(j,i)/=U.at(i,i);
        }
    }
    U.at(n,n)=A.at(n,n);
    for(unsigned i(0);i!=n;++i)
        U.at(n,n)-=L.at(n,i)*U.at(i,n);
    if(abs(U.at(n,n))<e)return true;
    return false;
}
}
E3t2::E3t2(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E3t2)
{
    ui->setupUi(this);
    ui->widget->setTitle("待分解矩阵");
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        if(ui->widget->x.empty()||ui->widget->x.n_cols!=ui->widget->x.n_rows)
        {
            QMessageBox::critical(parent,"错误","输入有误");
            return;
        }
        mat L,U;
        //这里需要lapack库, 若有lapack库配置可解开如下注释并修改CMake文件
//        if(!lu(L,U,ui->widget->x))
//        {
//            QMessageBox::critical(parent,"错误","输入有误");
//            return;
//        }
        if(Ex3::LU(L,U,ui->widget->x))
        {
            QMessageBox::critical(parent,"错误","输入有误");
            return;
        }
        ui->label_3->showMat(L);
        ui->label_4->showMat(U);
    });
#ifdef DEBUG32
    mat A{{6,2,1,-1},{2,4,1,0},{1,1,4,-1},{-1,0,-1,3}},L,U;
    Ex3::LU(L,U,A);
    L.print("L=");
    U.print("U=");
#endif
}

E3t2::~E3t2()
{
    delete ui;
}
