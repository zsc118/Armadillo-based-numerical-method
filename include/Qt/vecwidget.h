#ifndef VECWIDGET_H
#define VECWIDGET_H

#include <QWidget>
#include <armadillo>
namespace Ui {
class vecWidget;
}

class vecWidget : public QWidget
{
    Q_OBJECT

public:
    arma::vec x;
    explicit vecWidget(QWidget *parent = nullptr);
    ~vecWidget();
    void setTitle(const char*);

private:
    Ui::vecWidget *ui;
};

#endif // VECWIDGET_H
