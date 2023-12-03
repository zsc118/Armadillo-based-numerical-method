#ifndef E1T1_H
#define E1T1_H

#include <QWidget>

namespace Ui {
class E1t1;
}

class E1t1 : public QWidget
{
    Q_OBJECT

public:
    explicit E1t1(QWidget *parent = nullptr);
    ~E1t1();

private:
    Ui::E1t1 *ui;
};

#endif // E1T1_H
