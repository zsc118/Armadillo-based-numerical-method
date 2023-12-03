#ifndef SHOWWIDGET_H
#define SHOWWIDGET_H

#include <QDialog>
#include <armadillo>
namespace Ui {
class showWidget;
}

class showWidget : public QDialog
{
    Q_OBJECT

public:
    explicit showWidget(QWidget *parent,arma::mat& x);
    ~showWidget();

private:
    Ui::showWidget *ui;
};

#endif // SHOWWIDGET_H
