#ifndef MATLABEL_H
#define MATLABEL_H
#include <QLabel>
#include <armadillo>
// 用于显示矩阵
class matLabel:public QLabel
{
    Q_OBJECT
    arma::mat a;
public:
    matLabel(QWidget*);
    // 打印矩阵
    void showMat(const arma::mat&);
    void mousePressEvent(QMouseEvent*);
};

#endif // MATLABEL_H
