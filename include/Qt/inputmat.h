#ifndef INPUTMAT_H
#define INPUTMAT_H

#include <QDialog>
#include <armadillo>

namespace Ui {
class inputMat;
}

class inputMat : public QDialog
{
    Q_OBJECT
public:
    // m: 输入矩阵的行数
    // n: 输入矩阵的列数
    explicit inputMat(QWidget *parent, unsigned m, unsigned n, arma::mat& x);
    ~inputMat();

private:
    Ui::inputMat *ui;
};

#endif // INPUTMAT_H
