#include "e2t5.h"
#include "ui_e2t5.h"
#include "dxs.h"
#include <vector>
#include <string>
using namespace std;
extern double calStr(string s);
typedef dxs<double,unsigned> Dxs;
namespace Ex2 {
/*
 * Sturm序列法
 * p: 多项式
 * e: 精度
 * r: 返回孤立区间
 *
 * 不检查零次多项式, 默认次数大于等于1
 */
bool sturm(const Dxs& p,vector<Interval>& r,double e=1e-6)
{
    r.clear();
    double* q(p.a),A(*q);
    unsigned k(p.n);
    do
        if(*++q>A)
            A=*q;
    while(--k);
}
}
E2t5::E2t5(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E2t5)
{
    ui->setupUi(this);
}

E2t5::~E2t5()
{
    delete ui;
}
