#ifndef INPUTMATROWCOL_H
#define INPUTMATROWCOL_H

#include <QDialog>
#include <armadillo>
namespace Ui {
class inputMatRowCol;
}

class inputMatRowCol : public QDialog
{
    Q_OBJECT

public:
    explicit inputMatRowCol(QWidget *parent,arma::mat& x);
    ~inputMatRowCol();

private:
    Ui::inputMatRowCol *ui;
};

#endif // INPUTMATROWCOL_H
