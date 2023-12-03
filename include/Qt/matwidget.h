#ifndef MATWIDGET_H
#define MATWIDGET_H

#include <QWidget>
#include <armadillo>
namespace Ui {
class matWidget;
}

class matWidget : public QWidget
{
    Q_OBJECT

public:
    arma::mat x;
    explicit matWidget(QWidget *parent = nullptr);
    ~matWidget();
    void setTitle(const char*);

private:
    Ui::matWidget *ui;
};

#endif // MATWIDGET_H
