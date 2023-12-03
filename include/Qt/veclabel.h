#ifndef VECLABEL_H
#define VECLABEL_H
#include <QLabel>
#include <armadillo>
// 用于显示向量
class vecLabel:public QLabel
{
    Q_OBJECT
public:
    vecLabel(QWidget*);
    // 打印向量
    void showVec(const arma::vec&,const char*);
};

#endif // VECLABEL_H
