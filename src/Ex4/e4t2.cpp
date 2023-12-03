#include "e4t2.h"
#include "ui_e4t2.h"
#include <QMessageBox>
#include <QFile>
#include <QInputDialog>
#include <QFileDialog>
#include <QStandardPaths>
#include <string>
#include <vector>
#ifdef DEBUG42
#include <QDebug>
#endif
using namespace arma;
using namespace std;
namespace _MYFUNCTION {
extern bool check_unique(const vec&);
}
namespace Ex4 {
/*
 * Lagrange型Hermite插值法
 * x0:观测点
 * y0:观测数据(第i列为i阶导数)
 * m :每个节点有效的导数阶数
 * x1:插值点
 * y1:预测值
 */
//void Hermite_Lagrange(const vec& x0,const mat& y0,const vector<unsigned>& m,const vec& x1,vec& y1)
//{
//    if(x1.empty())
//    {
//        y1.clear();
//        return;
//    }
//    if(x0.empty()||y0.empty())throw "观测集为空！";
//    if(m.size()!=x0.n_elem)throw "导数阶数向量与观测点向量大小不一致！";
//    if(x0.n_elem!=y0.n_rows)throw "观测点与观测数据向量大小不一致！";
//    if(_MYFUNCTION::check_unique(x0))throw "观测集包含重复元素！";
//    y1.zeros(size(x1));
//    const unsigned end(-1),m_max(y0.n_cols-1);
//    unsigned*fact(new unsigned[y0.n_cols]),*p(fact);//阶乘
//    *p=1;
//    for(unsigned i(0);i!=m_max;++p)
//        *(p+1)=++i**p;
//    for(unsigned i(0);i!=m.size();++i)
//    {
//        unsigned k(m.at(i));
//        if(k>m_max)
//        {
//            delete[]fact;
//            throw "观测数据缺失！";
//        }
//        ++(p=fact+k);
//        do
//        {
//            vec t(size(x1),fill::value(y0.at(i,k)/ *--p));
//            unsigned j(end);
//            while(++j!=i)
//        }while(--k!=end);
//    }
//    delete[]fact;
//}

//差商表
class Difference_quotient_table
{
    double* x;//观测值
    double* y;//差商值
    unsigned *m;//当前导数阶数
    const mat* Y;//观测数据
    const vector<unsigned>* M;//导数阶数
    unsigned k;//差商数
    unsigned step;//观测点的步长
    unsigned *pre_m;//当前导数阶数的前一个位置
    unsigned fact;//阶乘
    double* x_out;//输出当前x
//    static constexpr double eps=1e-100;
public:
    Difference_quotient_table(const vec* x0,const mat* y0,const vector<unsigned>* M):Y(y0),step(0),fact(1)
    {
        unsigned i(k=(this->M=M)->size()),*m_p(m=new unsigned[k]+k);
        auto M_p(M->cend());
        do
        {
            if((*--m=*--M_p)>y0->n_cols)
            {
                delete[](m-=--i);
                throw "导数数据缺失！";
            }
            k+=*m;
        }
        while(--i);
        y=new double[k]+k;
        x=new double[k]+k;
        const double* y_p(&y0->at(M->size())),*x_p(&x0->at(M->size()));//armadillo库的矩阵是按列存储的
        do
        {
            *--x=*--x_p,*--y=*--y_p;
            if(i=*--m_p)do *--x=*x_p,*--y=*y_p;while(--i);
        }while(m_p!=m);
        pre_m=m-1;
        --(x_out=x);
    }
    ~Difference_quotient_table()noexcept
    {
        delete[] y;
        delete[] x;
        delete[] m;
    }
    //进行差商
    bool diff()
    {
        if(!--k)return false;//只有1个元素, 不需要计算
        unsigned d_k(++step);
        fact*=step;
        unsigned* m_p(pre_m),i(0);
        double* p(y),*q(x),*q1(q+step);
        do
            if(*q!=*q1)//abs(*p-=*(p+1))>eps
                (*p-=*(p+1))/=*q-*q1;
            else
            {
                while(!*++m_p);//找到第1个未完全计算的重节点差商
                unsigned j(*m_p);
                i+=--*m_p;
                q+=*m_p;
                q1+=*m_p;
                const double t(Y->at(m_p-m,d_k)/fact);
                *p=t;
                while(--j)*++p=t;
            }
        while(++q,++q1,++p,++i!=k);
        return true;
    }
    double& get_diff()noexcept
    {
        return *y;
    }
    double& get_x()noexcept
    {
        return *++x_out;
    }
};
/*
 * Hermite插值法
 * x0:观测点
 * y0:观测数据(第i列为i阶导数)
 * m :每个节点有效的导数阶数
 * x1:插值点
 * y1:预测值
 */
void Hermite(const vec& x0,const mat& y0,const vector<unsigned>& m,const vec& x1,vec& y1)
{
    if(x1.empty())
    {
        y1.clear();
        return;
    }
    if(x0.empty()||y0.empty())throw "观测集为空！";
    if(m.size()!=x0.n_elem)throw "导数阶数向量与观测点向量大小不一致！";
    if(x0.n_elem!=y0.n_rows)throw "观测点与观测数据向量大小不一致！";
    if(_MYFUNCTION::check_unique(x0))throw "观测集包含重复元素！";
//    double*X(new double[y0.n_elem]),*p(X),*Y(new double[y0.n_elem]),*u(Y);//这里出于运行效率考虑, 可能多申请了一些空间, 但避免了一次求和操作
//    auto q(x0.cbegin());
//    auto m_p(m.cbegin());
//    auto v(y0.cbegin());
//    *p=*q,*u=*v;
//    const unsigned& m_0(*m_p);
//    unsigned n(m_0);
//    for(unsigned i(0);i!=n;++i)
//        *++p=*q,*++u=*v;
//    while(++q!=x0.cend())
//    {
//        ++v;
//        n+=*++m_p;
//        for(unsigned i(_MYFUNCTION::max_u);i!=*m_p;++i)
//            *++p=*q,*++u=*v;
//    }
//    y1.set_size(size(x1)).fill(*Y);
//    vec t(ones(size(x1)));
//    n+=m.size();
//    unsigned i(0),k(n);
//    do
//    {
//        t%=x1-X[i];
//        p=X,u=Y;
//        *u=i++<=m_0?y0.at(0,i):(*(u+1)-*u)/(*(p+i)-*p);
//        unsigned t_k(--k);
//    }while(k);
    Difference_quotient_table D(&x0,&y0,&m);
    y1.set_size(x1.n_elem).fill(D.get_diff());
    vec x(ones(x1.n_elem));
    while(D.diff())
        y1+=D.get_diff()*(x%=x1-D.get_x());
}
}
E4t2::E4t2(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E4t2)
{
    ui->setupUi(this);
#ifdef DEBUG42
    auto ans=[](double x)->double{
        return (((1.5*x-2)*x+1.5)*x+2)*x+1;
    };
    vec x0{0.,1.},x1(linspace(0,1)),f_x1(x1),y1;
    f_x1.transform(ans);
    mat y0{{1,2,3},{4,5,0}};
    vector<unsigned> m;
    m.push_back(2);
    m.push_back(1);
    try {
        Ex4::Hermite(x0,y0,m,x1,y1);
//        for(auto& i:y1)qDebug()<<i;
        QFile file("D:/1.csv");
        file.open(QIODevice::WriteOnly);
        for(unsigned i(0);i!=101;++i)
            file.write((to_string(x1.at(i))+','+to_string(y1.at(i))+','+to_string(f_x1.at(i))+'\n').c_str());
//        qDebug();
//        for(auto& i:f_x1)qDebug()<<i;
    } catch (const char* s) {
        qDebug()<<s;
    }
#else
    ui->widget->setTitle("观测点");
    ui->widget_2->setTitle("观测数据");
    ui->widget_3->setTitle("插值点");
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        if(ui->widget_3->x.empty())
        {
            QMessageBox::warning(this,"Hermite插值","插值集为空！");
            return;
        }
        vec t;
        try {
            vector<unsigned> m;
            unsigned i(0);
            while(i!=ui->widget->x.n_elem)
            {
                string s("请输入第");
                int t(QInputDialog::getInt(this,"Hermite插值",((s+=to_string(++i))+="个未知数所含的导数的最高阶数:").c_str()));
                if(t<0)throw "导数阶数输入有误！";
                m.push_back(unsigned(t));
            }
            Ex4::Hermite(ui->widget->x,ui->widget_2->x,m,ui->widget_3->x,t);
        } catch (const char* s) {
            QMessageBox::critical(this,"Hermite插值",s);
            return;
        }
        QFile file(QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"));
        if(!file.open(QIODevice::WriteOnly))
        {
            QMessageBox::critical(this,"Hermite插值","文件保存失败！");
            return;
        }
        file.write("x,y");
        for(unsigned i(0);i!=t.n_elem;++i)
        {
            string s("\n");
            file.write((((s+=to_string(ui->widget_3->x.at(i)))+=',')+=to_string(t.at(i))).c_str());
        }
        file.close();
        QMessageBox::information(this,"Hermite插值","导出成功！");
    });
#endif
}

E4t2::~E4t2()
{
    delete ui;
}
