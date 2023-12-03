#ifndef E4T3_1_H
#define E4T3_1_H
#include <armadillo>
#include <QDialog>
namespace Ex4 {
//边界条件
struct boundary_conditions
{
    unsigned char t;//边界条件类型
    //第1位: 0:m/M型 1:周期边界条件
    //第2位: 0:第0个节点取1阶导 1:第0个节点取2阶导
    //第3位: 0:第n个节点取1阶导 1:第n个节点取2阶导
    double y1,y2;
    //默认构造:周期边界
    boundary_conditions():t('\1'){}
    boundary_conditions(unsigned char k,const double& w1=0.0,const double& w2=0.0):t(k),y1(w1),y2(w2){}
};
//Lagrange型边界条件
//不检查x,y的维数匹配, x元素是否唯一及元素个数问题
//公式推导请参考Lagrange_boundary.ggb
boundary_conditions& Lagrange_boundary(const arma::vec&,const arma::mat&,boundary_conditions&)noexcept;
}
namespace Ui {
class E4t3_1;
}

class E4t3_1 : public QDialog
{
    Q_OBJECT

public:
    explicit E4t3_1(QWidget *parent,Ex4::boundary_conditions&,const arma::vec&,const arma::mat&);
    ~E4t3_1();

private:
    Ui::E4t3_1 *ui;
};

#endif // E4T3_1_H
