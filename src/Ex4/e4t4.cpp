#include "e4t4.h"
#include "ui_e4t4.h"
#include <stdio.h>
#include <QMessageBox>
#include <QFileDialog>
#include <QStandardPaths>
using namespace arma;
#ifdef DEBUG44
#include <QDebug>
#endif
namespace _MYFUNCTION {
#ifdef DEBUG44
//产生2维均匀网格
void generate_2grid(mat& X,double xa,double xb,double ya,double yb,unsigned n=100)
{
    if(xa>xb)
    {
        double t(xa);
        xa=xb;
        xb=t;
    }
    if(ya>yb)
    {
        double t(ya);
        ya=yb;
        yb=t;
    }
    unsigned n0(n),i(++n);
    X.set_size(2,n*n);
    double* p(&X.at(0)),xstep((xb-xa)/n0),ystep((yb-ya)/n0),x(xa),y(ya);
    --p;
    do
    {
        unsigned j(n);
        do
        {
            *++p=x;
            *++p=y;
            y+=ystep;
        }while(--j);
        x+=xstep;
        y=ya;
    }while(--i);
}
#endif
/*
 * 检查矩阵是否包含重复列向量
 *  包含  :true
 *  不包含:false
 */
bool check_unique(const mat& x)
{
    unsigned i(x.n_cols);
    do
    {
        unsigned j(--i);
        while(j)if(all(x.col(i)==x.col(--j)))return true;
    }while(i);
    return false;
}
}
namespace Ex4 {
/*
 * 最邻近内插法
 * x0:观测点
 * y0:观测数据
 * x1:插值点
 * y1:预测值
 */
void nearest_neighbor(const mat& x0,const mat& y0,const mat& x1,mat& y1)
{
    if(x1.empty())
    {
        y1.clear();
        return;
    }
    if(x0.empty()||y0.empty())throw "观测集为空！";
    if(x0.n_cols!=y0.n_cols)throw "观测点与观测数据向量大小不一致！";
    if(x1.n_rows!=x0.n_rows)throw "插值点与观测点维数不同！";
    if(x0.n_cols==1)
    {
        y1=repmat(y0,1,x1.n_cols);
        return;
    }
    if(_MYFUNCTION::check_unique(x0))throw "观测集包含重复元素！";
    unsigned i(x1.n_cols);
    y1.set_size(y0.n_rows,i);
    do
    {
        auto x(x1.col(--i));
        unsigned j(x0.n_cols),k(--j);
        double m(norm(x-x0.col(j)));
        while(j)
        {
            double t(norm(x-x0.col(--j)));
            if(t<m)
                m=t,k=j;
        }
        y1.col(i)=y0.col(k);
    }while(i);
}
/*
 * 基于TIN的自然邻居法
 * x0:观测点
 * y0:观测数据
 * x1:插值点
 * y1:预测值
 */
void TIN_neighbor(const mat& x0,const mat& y0,const mat& x1,mat& y1)
{
    if(x1.empty())
    {
        y1.clear();
        return;
    }
    if(x0.empty()||y0.empty())throw "观测集为空！";
    if(x0.n_cols!=y0.n_cols)throw "观测点与观测数据向量大小不一致！";
    if(x1.n_rows!=x0.n_rows)throw "插值点与观测点维数不同！";
    if(x0.n_rows>=x0.n_cols)throw "观测点数量过少！";
    if(_MYFUNCTION::check_unique(x0))throw "观测集包含重复元素！";
    unsigned i(x1.n_cols),n1(x0.n_rows),nd(n1++),n(nd--);
    y1.set_size(y0.n_rows,i);
    mat R(n,n1,fill::none);
    unsigned j0(x0.n_cols);
    vec dis(j0,fill::none);
    double* p0(&dis.at(--j0));
    do
    {
        auto x(x1.col(--i)),y(y1.col(i));
        unsigned j(j0);
        double* p(p0);
        *p=norm(x0.col(j)-x);
        while(j)*--p=norm(x0.col(--j)-x);
        uvec index(sort_index(dis));
        auto x_s0(x0.col(index.at(0)));
        R.col(0)=x-x_s0;
        R.col(j=1)=x0.col(index.at(1))-x_s0;
        while(rank(R.cols(0,j))>j)
        {
            ++j;
            R.col(j)=x0.col(index.at(j))-x_s0;
        }
        vec K(solve(R.cols(1,j),R.col(0)));
        auto y_s0(y0.col(index.at(0)));
        y=y_s0;
        while(j--)y+=K.at(j)*(y0.col(index.at(j))-y_s0);
    }while(i);
}
/*
 * 反距离加权插值法
 * x0:观测点
 * y0:观测数据
 * x1:插值点
 * y1:预测值
 */
void inverse_distance_weighting(const mat& x0,const mat& y0,const mat& x1,mat& y1)
{
    if(x1.empty())
    {
        y1.clear();
        return;
    }
    if(x0.empty()||y0.empty())throw "观测集为空！";
    if(x0.n_cols!=y0.n_cols)throw "观测点与观测数据向量大小不一致！";
    if(x1.n_rows!=x0.n_rows)throw "插值点与观测点维数不同！";
    if(x0.n_cols==1)
    {
        y1=repmat(y0,1,x1.n_cols);
        return;
    }
    if(_MYFUNCTION::check_unique(x0))throw "观测集包含重复元素！";
    unsigned i(x1.n_cols),n(x0.n_cols);
    y1.set_size(y0.n_rows,i);
    double *a(new double[n]),*p0(a+n);
    do
    {
        auto x(x1.col(--i)),y(y1.col(i));
        unsigned j(n);
        double* p(p0),sum(0);
        bool flag(false);
        do
        {
            if(!(*--p=norm(x0.col(--j)-x)))
            {
                flag=true;
                break;
            }
            sum+=*p=1/ *p;
        }while(j);
        if(flag)
        {
            y1.col(i)=y0.col(j);
            continue;
        }
        *p/=sum;
        while(++p!=p0)*p/=sum;
        y=*p*y0.col(--(j=n));
        while(j)y+=*--p*y0.col(--j);
    }while(i);
    delete[]a;
}
}
E4t4::E4t4(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E4t4)
{
    ui->setupUi(this);
#ifdef DEBUG44
    mat x0{{2,2,9,9},{1,6,1,6}},y0{60,55,57.5,70},x1,y1;
    _MYFUNCTION::generate_2grid(x1,0,10,0,10);
    try {
//        Ex4::nearest_neighbor(x0,y0,x1,y1);
//        Ex4::TIN_neighbor(x0,y0,x1,y1);
        Ex4::inverse_distance_weighting(x0,y0,x1,y1);
        FILE* file;
        if(fopen_s(&file,"D:/1.csv","w"))throw "文件保存失败！";
        for(unsigned i(0);i!=y1.n_cols;++i)
            fprintf(file,"%.14f,%.14f,%.14f\n",x1.at(0,i),x1.at(1,i),y1.at(i));
    } catch (const char* s) {
        qDebug()<<s;
    }
#else
    ui->widget->setTitle("观测点");
    ui->widget_2->setTitle("观测数据");
    ui->widget_3->setTitle("插值点");
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        mat &x0(ui->widget->x),&y0(ui->widget_2->x),&x1(ui->widget_3->x);
        if(x1.empty())
        {
            QMessageBox::warning(this,"高维插值","插值点为空！");
            return;
        }
        try {
            mat y1;
            FILE* file;
            if(fopen_s(&file,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()").toLatin1().data(),"w"))
                throw "文件保存失败！";
            if(ui->radioButton->isChecked())
                Ex4::nearest_neighbor(x0,y0,x1,y1);
            else if(ui->radioButton_2->isChecked())
                Ex4::TIN_neighbor(x0,y0,x1,y1);
            else
                Ex4::inverse_distance_weighting(x0,y0,x1,y1);
            unsigned n(y1.n_rows-1);
            for(unsigned i(0);i!=x1.n_cols;++i)
            {
                for(unsigned j(0);j!=x1.n_rows;++j)
                    fprintf(file,"%.14f,",x1.at(j,i));
                for(unsigned j(0);j!=n;++j)
                    fprintf(file,"%.14f,",y1.at(j,i));
                fprintf(file,"%.14f\n",y1.at(n,i));
            }
            fclose(file);
        } catch (const char* s) {
            QMessageBox::critical(this,"高维插值",s);
            return;
        }
    });
#endif
}

E4t4::~E4t4()
{
    delete ui;
}
