#ifndef E3T1_H
#define E3T1_H

#include <QWidget>

namespace Ui {
class E3t1;
}

class E3t1 : public QWidget
{
    Q_OBJECT

public:
    explicit E3t1(QWidget *parent = nullptr);
    ~E3t1();

private:
    Ui::E3t1 *ui;
};

#endif // E3T1_H
