#include "e3t1.h"
#include "ui_e3t1.h"
#include <armadillo>
#include <vector>
#include <QFile>
#include <QStandardPaths>
#include <QMessageBox>
#include <string>
#include <QFileDialog>
#if defined(PRINT_LATEX)||defined(DEBUG31)
#include <QDebug>
#endif
#ifdef ZSC_TIME
#include <time.h>
#include <QDebug>
#endif
using namespace arma;
using namespace std;
namespace _MYFUNCTION {
#ifdef PRINT_LATEX
void print_latex_mat1(const mat& A)
{
    qDebug()<<"\\begin{pmatrix}";
    for(unsigned i(0);i<A.n_rows;++i)
    {
        string s;
        for(unsigned j(0);j<A.n_rows;++j)
            (s+=to_string(A.at(i,j)))+='&';
        qDebug()<<s.c_str()<<A.at(i,A.n_rows)<<"\\\\";
    }
    qDebug()<<"\\end{pmatrix}\\\\\\rightarrow&";
}
#endif
void print_mat1(QFile& file,const mat& A)
{
    string s;
    for(unsigned i(0);i<A.n_rows;++i)
    {
        for(unsigned j(0);j<A.n_rows;++j)
            (s+=to_string(A.at(i,j)))+=',';
        (s+=to_string(A.at(i,A.n_rows)))+='\n';
    }
    file.write((s+='\n').c_str());
}
}
namespace Ex3 {
/*
 * 顺序Gauss消去法
 * A : 系数矩阵
 * b : 常数向量
 * uA: 系数矩阵消元过程
 * ub: 常数向量消元过程
 * e : 判零标准
 *
 * 返回(bool):
 *  true : 求解失败
 *  false: 求解成功
 *
 * 不检查矩阵维数是否匹配问题
 */
bool Sequential_Gauss(mat A,vec b,vector<mat>& uA,vector<vec>& ub,const double& e=1e-6)
{
    uA.clear();
    ub.clear();
    uA.push_back(A);
    ub.push_back(b);
    unsigned n(A.n_cols-1),m(-1);
    //消元
    for(unsigned i(0);i!=n;++i)//arma::mat的迭代器不太好用, 最后还是用了下标形式, 效率相比迭代器低一些
    {
        if(A.at(i,i)<e&&A.at(i,i)>-e)return true;
        for(unsigned j(i+1);j!=A.n_cols;++j)
        {
            double t(A.at(j,i)/A.at(i,i));
            A.at(j,i)=0;
            b.at(j)-=t*b.at(i);
            for(unsigned k(i+2);k<A.n_cols;++k)
                A.at(j,k)-=t*A.at(i,k);
        }
        uA.push_back(A);
        ub.push_back(b);
    }
    //回代
    do
    {
        if(A.at(n,n)<e&&A.at(n,n)>-e)return true;
        b.at(n)/=A.at(n,n);
        A.at(n,n)=1;
        for(unsigned i(0);i!=n;++i)
        {
            b.at(i)-=A.at(i,n)*b.at(n);
            A.at(i,n)=0;
        }
        uA.push_back(A);
        ub.push_back(b);
    }while(--n!=m);
    return false;
}
/*
 * 顺序Gauss消去法
 * A   : 增广矩阵
 * path: 迭代过程保存
 * e   : 判零标准
 * res : 结果
 *
 * 返回(bool):
 *  true : 求解失败
 *  false: 求解成功
 */
bool Sequential_Gauss(mat A,const QString& path,const double& e,vec& res)
{
    if(A.n_cols!=A.n_rows+1)return true;
    QFile file(path);
    if(!file.open(QIODevice::WriteOnly))return true;
    _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
    _MYFUNCTION::print_latex_mat1(A);
#endif
    unsigned n(A.n_rows-1),m(-1);
    for(unsigned i(0);i!=n;++i)//消元
    {
        if(A.at(i,i)<e&&A.at(i,i)>-e)return true;
        for(unsigned j(i+1);j!=A.n_rows;++j)
        {
            double t(A.at(j,i)/A.at(i,i));
            A.at(j,i)=0;
            for(unsigned k(i+1);k<A.n_cols;++k)
                A.at(j,k)-=t*A.at(i,k);
        }
        _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
        _MYFUNCTION::print_latex_mat1(A);
#endif
    }
    res.zeros(A.n_rows,1);
    do//回代
    {
        if(A.at(n,n)<e&&A.at(n,n)>-e)return true;
        res.at(n)=A.at(n,A.n_rows)/=A.at(n,n);
        A.at(n,n)=1;
        for(unsigned i(0);i!=n;++i)
        {
            A.at(i,A.n_rows)-=A.at(i,n)*A.at(n,A.n_rows);
            A.at(i,n)=0;
        }
        _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
        _MYFUNCTION::print_latex_mat1(A);
#endif
    }while(--n!=m);
    return false;
}
/*
 * 列选主元Gauss消去法
 * A : 系数矩阵
 * b : 常数向量
 * uA: 系数矩阵消元过程
 * ub: 常数向量消元过程
 * e : 判零标准
 *
 * 返回(bool):
 *  true : 求解失败
 *  false: 求解成功
 *
 * 不检查矩阵维数是否匹配问题
 */
bool Col_Gauss(mat A,vec b,vector<mat>& uA,vector<vec>& ub,const double& e=1e-6)
{
    uA.clear();
    ub.clear();
    uA.push_back(A);
    ub.push_back(b);
    unsigned n(A.n_cols-1),m(-1);
    for(unsigned i(0);i!=n;++i)//消元
    {
        unsigned max(i);//选列主元
        for(unsigned j(i+1);j!=A.n_cols;++j)
            if(abs(A.at(j,i))>abs(A.at(max,i)))
                max=j;
        if(abs(A.at(max,i))<e)return true;
        if(max!=i)
        {
            for(unsigned j(i);j!=A.n_cols;++j)
            {
                double t(A.at(i,j));
                A.at(i,j)=A.at(max,j);
                A.at(max,j)=t;
            }
            double t(b.at(max));
            b.at(max)=b.at(i);
            b.at(i)=t;
            uA.push_back(A);
            ub.push_back(b);
        }
        for(unsigned j(i+1);j!=A.n_cols;++j)
        {
            double t(A.at(j,i)/A.at(i,i));
            A.at(j,i)=0;
            b.at(j)-=t*b.at(i);
            for(unsigned k(i+1);k<A.n_cols;++k)
                A.at(j,k)-=t*A.at(i,k);
        }
        uA.push_back(A);
        ub.push_back(b);
    }
    do//回代
    {
        if(A.at(n,n)<e&&A.at(n,n)>-e)return true;
        b.at(n)/=A.at(n,n);
        A.at(n,n)=1;
        for(unsigned i(0);i!=n;++i)
        {
            b.at(i)-=A.at(i,n)*b.at(n);
            A.at(i,n)=0;
        }
        uA.push_back(A);
        ub.push_back(b);
    }while(--n!=m);
    return false;
}
/*
 * 列选主元Gauss消去法
 * A   : 增广矩阵
 * path: 迭代过程保存
 * e   : 判零标准
 * res : 结果
 *
 * 返回(bool):
 *  true : 求解失败
 *  false: 求解成功
 */
bool Col_Gauss(mat A,const QString& path,const double& e,vec& res)
{
    if(A.n_cols!=A.n_rows+1)return true;
    QFile file(path);
    if(!file.open(QIODevice::WriteOnly))return true;
    _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
    _MYFUNCTION::print_latex_mat1(A);
#endif
    unsigned n(A.n_rows-1),m(-1);
    for(unsigned i(0);i!=n;++i)//消元
    {
        unsigned max(i);//选列主元
        for(unsigned j(i+1);j!=A.n_rows;++j)
            if(abs(A.at(j,i))>abs(A.at(max,i)))
                max=j;
        if(abs(A.at(max,i))<e)return true;
        if(max!=i)
        {
            for(unsigned j(i);j!=A.n_cols;++j)
            {
                double t(A.at(i,j));
                A.at(i,j)=A.at(max,j);
                A.at(max,j)=t;
            }
            _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
            _MYFUNCTION::print_latex_mat1(A);
#endif
        }
        for(unsigned j(i+1);j!=A.n_rows;++j)
        {
            double t(A.at(j,i)/A.at(i,i));
            A.at(j,i)=0;
            for(unsigned k(i+1);k<A.n_cols;++k)
                A.at(j,k)-=t*A.at(i,k);
        }
        _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
        _MYFUNCTION::print_latex_mat1(A);
#endif
    }
    res.zeros(A.n_rows,1);
    do//回代
    {
        if(A.at(n,n)<e&&A.at(n,n)>-e)return true;
        res.at(n)=A.at(n,A.n_rows)/=A.at(n,n);
        A.at(n,n)=1;
        for(unsigned i(0);i!=n;++i)
        {
            A.at(i,A.n_rows)-=A.at(i,n)*A.at(n,A.n_rows);
            A.at(i,n)=0;
        }
        _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
        _MYFUNCTION::print_latex_mat1(A);
#endif
    }while(--n!=m);
    return false;
}
/*
 * 全选主元Gauss消去法
 * A : 系数矩阵
 * b : 常数向量
 * uA: 系数矩阵消元过程
 * ub: 常数向量消元过程
 * r : 解向量
 * e : 判零标准
 *
 * 返回(bool):
 *  true : 求解失败
 *  false: 求解成功
 *
 * 不检查矩阵维数是否匹配问题
 */
bool All_Gauss(mat A,vec b,vector<mat>& uA,vector<vec>& ub,vec& r,const double& e=1e-6)
{
    uA.clear();
    ub.clear();
    uA.push_back(A);
    ub.push_back(b);
    r.zeros(A.n_cols);
    unsigned n(A.n_cols-1),m(-1),*a(new unsigned[A.n_cols]);
    for(unsigned t(0);t<A.n_rows;++t)//初始化排序向量
        a[t]=t;
    for(unsigned i(0);i!=n;++i)//消元
    {
        unsigned x(i),y(i);//选主元
        for(unsigned j(i);j!=A.n_cols;++j)
            for(unsigned k(i);k!=A.n_cols;++k)
                if(abs(A.at(j,k))>abs(A.at(x,y)))
                {
                    x=j;
                    y=k;
                }
        if(abs(A.at(x,y))<e)
        {
            delete[]a;
            return true;
        }
        bool flag(false);
        if(x!=i)
        {
            for(unsigned j(i);j<A.n_cols;++j)
            {
                double t(A.at(x,j));
                A.at(x,j)=A.at(i,j);
                A.at(i,j)=t;
            }
            double t(b.at(x));
            b.at(x)=b.at(i);
            b.at(i)=t;
            flag=true;
        }
        if(y!=i)
        {
            for(unsigned j(0);j<A.n_cols;++j)
            {
                double t(A.at(j,i));
                A.at(j,i)=A.at(j,y);
                A.at(j,y)=t;
            }
            flag=true;
            unsigned t(a[i]);
            a[i]=a[y];
            a[y]=t;
        }
        if(flag)
        {
            uA.push_back(A);
            ub.push_back(b);
        }
        for(unsigned j(i+1);j!=A.n_cols;++j)
        {
            double t(A.at(j,i)/A.at(i,i));
            A.at(j,i)=0;
            b.at(j)-=t*b.at(i);
            for(unsigned k(i+2);k<A.n_cols;++k)
                A.at(j,k)-=t*A.at(i,k);
        }
        uA.push_back(A);
        ub.push_back(b);
    }
    do//回代
    {
        if(abs(A.at(n,n))<e)
        {
            delete[]a;
            return true;
        }
        r.at(a[n])=b.at(n)/=A.at(n,n);
        A.at(n,n)=1;
        for(unsigned i(0);i!=n;++i)
        {
            b.at(i)-=A.at(i,n)*b.at(n);
            A.at(i,n)=0;
        }
        uA.push_back(A);
        ub.push_back(b);
    }while(--n!=m);
    delete[]a;
    return false;
}
/*
 * 全选主元Gauss消去法
 * A   : 增广矩阵
 * path: 迭代过程保存
 * e   : 判零标准
 * res : 结果
 *
 * 返回(bool):
 *  true : 求解失败
 *  false: 求解成功
 */
bool All_Gauss(mat A,const QString& path,const double& e,vec& res)
{
    if(A.n_cols!=A.n_rows+1)return true;
    QFile file(path);
    if(!file.open(QIODevice::WriteOnly))return true;
    _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
    _MYFUNCTION::print_latex_mat1(A);
#endif
    unsigned n(A.n_rows-1),m(-1),*a(new unsigned[A.n_rows]);
    for(unsigned t(0);t<A.n_rows;++t)//初始化排序向量
        a[t]=t;
    for(unsigned i(0);i!=n;++i)//消元
    {
        unsigned x(i),y(i);//选主元
        for(unsigned j(i);j!=A.n_cols;++j)
            for(unsigned k(i);k!=A.n_cols;++k)
                if(abs(A.at(j,k))>abs(A.at(x,y)))
                {
                    x=j;
                    y=k;
                }
        if(abs(A.at(x,y))<e)
        {
            delete[]a;
            return true;
        }
        bool flag(false);
        if(x!=i)
        {
            for(unsigned j(i);j!=A.n_cols;++j)
            {
                double t(A.at(x,j));
                A.at(x,j)=A.at(i,j);
                A.at(i,j)=t;
            }
            flag=true;
        }
        if(y!=i)
        {
            for(unsigned j(0);j!=A.n_rows;++j)
            {
                double t(A.at(j,i));
                A.at(j,i)=A.at(j,y);
                A.at(j,y)=t;
            }
            flag=true;
            unsigned t(a[i]);
            a[i]=a[y];
            a[y]=t;
        }
        if(flag)
        {
            _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
            _MYFUNCTION::print_latex_mat1(A);
#endif
        }
        for(unsigned j(i+1);j!=A.n_rows;++j)
        {
            double t(A.at(j,i)/A.at(i,i));
            A.at(j,i)=0;
            for(unsigned k(i+1);k<A.n_cols;++k)
                A.at(j,k)-=t*A.at(i,k);
        }
        _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
        _MYFUNCTION::print_latex_mat1(A);
#endif
    }
    res.zeros(A.n_rows,1);
    do//回代
    {
        if(abs(A.at(n,n))<e)
        {
            delete[]a;
            return true;
        }
        res.at(a[n])=A.at(n,A.n_rows)/=A.at(n,n);
        A.at(n,n)=1;
        for(unsigned i(0);i!=n;++i)
        {
            A.at(i,A.n_rows)-=A.at(i,n)*A.at(n,A.n_rows);
            A.at(i,n)=0;
        }
        _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
        _MYFUNCTION::print_latex_mat1(A);
#endif
    }while(--n!=m);
    delete[]a;
    return false;
}
/*
 * 列选主元Gauss-Jordan消去法
 * A : 系数矩阵
 * b : 常数向量
 * uA: 系数矩阵消元过程
 * ub: 常数向量消元过程
 * e : 判零标准
 *
 * 返回(bool):
 *  true : 求解失败
 *  false: 求解成功
 *
 * 不检查矩阵维数是否匹配问题
 */
bool Gauss_Jordan(mat A,vec b,vector<mat>& uA,vector<vec>& ub,const double& e=1e-6)
{
    uA.clear();
    ub.clear();
    uA.push_back(A);
    ub.push_back(b);
    unsigned n(A.n_cols-1);
    for(unsigned i(0);i!=n;++i)
    {
        unsigned max(i);//选列主元
        for(unsigned j(i+1);j!=A.n_cols;++j)
            if(abs(A.at(j,i))>abs(A.at(max,i)))
                max=j;
        if(abs(A.at(max,i))<e)return true;
        if(max!=i)
        {
            for(unsigned j(i);j!=A.n_cols;++j)
            {
                double t(A.at(i,j));
                A.at(i,j)=A.at(max,j);
                A.at(max,j)=t;
            }
            double t(b.at(max));
            b.at(max)=b.at(i);
            b.at(i)=t;
            uA.push_back(A);
            ub.push_back(b);
        }
        b.at(i)/=A.at(i,i);
        for(unsigned k(i+1);k!=A.n_cols;++k)
            A.at(i,k)/=A.at(i,i);
        for(unsigned j(i+1);j<A.n_rows;++j)
        {
            b.at(j)-=b.at(i)*A.at(j,i);
            for(unsigned k(i+1);k<A.n_cols;++k)
                A.at(j,k)-=A.at(i,k)*A.at(j,i);
            A.at(j,i)=0;
        }
        for(unsigned j(0);j!=i;++j)
        {
            b.at(j)-=b.at(i)*A.at(j,i);
            for(unsigned k(i+1);k<A.n_cols;++k)
                A.at(j,k)-=A.at(i,k)*A.at(j,i);
            A.at(j,i)=0;
        }
        A.at(i,i)=1;
        uA.push_back(A);
        ub.push_back(b);
    }
    return false;
}
/*
 * 列选主元Gauss-Jordan消去法
 * A   : 增广矩阵
 * path: 迭代过程保存
 * e   : 判零标准
 * res : 结果
 *
 * 返回(bool):
 *  true : 求解失败
 *  false: 求解成功
 */
bool Gauss_Jordan(mat A,const QString& path,const double& e,vec& res)
{
    if(A.n_cols!=A.n_rows+1)return true;
    QFile file(path);
    if(!file.open(QIODevice::WriteOnly))return true;
    _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
    _MYFUNCTION::print_latex_mat1(A);
#endif
    res.zeros(A.n_rows);
    for(unsigned i(0);i!=A.n_rows;++i)
    {
        unsigned max(i);//选列主元
        for(unsigned j(i+1);j!=A.n_rows;++j)
            if(abs(A.at(j,i))>abs(A.at(max,i)))
                max=j;
        if(abs(A.at(max,i))<e)return true;
        if(max!=i)
        {
            for(unsigned j(i);j!=A.n_cols;++j)
            {
                double t(A.at(i,j));
                A.at(i,j)=A.at(max,j);
                A.at(max,j)=t;
            }
            _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
            _MYFUNCTION::print_latex_mat1(A);
#endif
        }
        for(unsigned k(i+1);k!=A.n_cols;++k)
            A.at(i,k)/=A.at(i,i);
        for(unsigned j(i+1);j<A.n_rows;++j)
        {
            for(unsigned k(i+1);k!=A.n_cols;++k)
                A.at(j,k)-=A.at(i,k)*A.at(j,i);
            A.at(j,i)=0;
        }
        for(unsigned j(0);j!=i;++j)
        {
            for(unsigned k(i+1);k<A.n_cols;++k)
                A.at(j,k)-=A.at(i,k)*A.at(j,i);
            A.at(j,i)=0;
        }
        A.at(i,i)=1;
        _MYFUNCTION::print_mat1(file,A);
#ifdef PRINT_LATEX
        _MYFUNCTION::print_latex_mat1(A);
#endif
    }
    return false;
}
}
E3t1::E3t1(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E3t1)
{
#ifdef DEBUG31
//    mat A{{0.101,2.304,3.555,1.183},{-1.347,3.712,4.623,2.137},{-2.835,1.072,5.643,3.035}};
    mat A{{0.101,2.304,3.555},{-1.347,3.712,4.623},{-2.835,1.072,5.643}};
    vec b{1.183,2.137,3.035};
//    vec t;
//    if(Ex3::Col_Gauss(A,"D:/1.csv",1e-6,t))qDebug()<<0;
    vec ans(solve(A,b));
    ans.print("ans=");
#endif
    ui->setupUi(this);
    ui->widget->setTitle("系数矩阵");
    ui->widget_2->setTitle("常数向量");
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        if(ui->widget->x.empty()||ui->widget_2->x.empty())
        {
            QMessageBox::critical(parent,"错误","输入为空！");
            return;
        }
        vec r;
        if(ui->radioButton->isChecked())
            if(Ex3::Sequential_Gauss(join_rows(ui->widget->x,ui->widget_2->x),QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"),1e-6,r))
            {
                QMessageBox::critical(parent,"错误","输入有误！");
                r.clear();
            }
            else
                QMessageBox::information(parent,"Gauss消元法","导出成功！");
        else if(ui->radioButton_2->isChecked())
            if(Ex3::Col_Gauss(join_rows(ui->widget->x,ui->widget_2->x),QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"),1e-6,r))
            {
                QMessageBox::critical(parent,"错误","输入有误！");
                r.clear();
            }
            else
            {
                QMessageBox::information(parent,"Gauss消元法","导出成功！");
            }
        else if(ui->radioButton_3->isChecked())
            if(Ex3::All_Gauss(join_rows(ui->widget->x,ui->widget_2->x),QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"),1e-6,r))
            {
                QMessageBox::critical(parent,"错误","输入有误！");
                r.clear();
            }
            else
                QMessageBox::information(parent,"Gauss消元法","导出成功！");
        else
            if(Ex3::Gauss_Jordan(join_rows(ui->widget->x,ui->widget_2->x),QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"),1e-6,r))
            {
                QMessageBox::critical(parent,"错误","输入有误！");
                r.clear();
            }
            else
                QMessageBox::information(parent,"Gauss消元法","导出成功！");
        ui->label->showVec(r,"解向量");
#ifdef ZSC_TIME
        unsigned short end(-1);
        {
            clock_t start_time=clock();
            for(unsigned short i(0);i!=end;++i)Ex3::Sequential_Gauss(join_rows(ui->widget->x,ui->widget_2->x),"D:/1.csv",1e-6,r);
            clock_t end_time=clock();
            qDebug()<<"顺序Gauss:"<<(end_time-start_time)/65536.0;
        }
        {
            clock_t start_time=clock();
            for(unsigned short i(0);i!=end;++i)Ex3::Col_Gauss(join_rows(ui->widget->x,ui->widget_2->x),"D:/1.csv",1e-6,r);
            clock_t end_time=clock();
            qDebug()<<"列主元Gauss:"<<(end_time-start_time)/65536.0;
        }
        {
            clock_t start_time=clock();
            for(unsigned short i(0);i!=end;++i)Ex3::All_Gauss(join_rows(ui->widget->x,ui->widget_2->x),"D:/1.csv",1e-6,r);
            clock_t end_time=clock();
            qDebug()<<"全主元Gauss:"<<(end_time-start_time)/65536.0;
        }
        {
            clock_t start_time=clock();
            for(unsigned short i(0);i!=end;++i)Ex3::Gauss_Jordan(join_rows(ui->widget->x,ui->widget_2->x),"D:/1.csv",1e-6,r);
            clock_t end_time=clock();
            qDebug()<<"Gauss-Jordan:"<<(end_time-start_time)/65536.0;
        }
#endif
    });
}

E3t1::~E3t1()
{
    delete ui;
}
