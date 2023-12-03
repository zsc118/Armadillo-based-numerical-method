#ifndef E4T1_H
#define E4T1_H

#include <QWidget>

namespace Ui {
class E4t1;
}

class E4t1 : public QWidget
{
    Q_OBJECT

public:
    explicit E4t1(QWidget *parent = nullptr);
    ~E4t1();

private:
    Ui::E4t1 *ui;
};

#endif // E4T1_H
