#ifndef E6T1_H
#define E6T1_H

#include <QWidget>

namespace Ui {
class E6t1;
}

class E6t1 : public QWidget
{
    Q_OBJECT

public:
    explicit E6t1(QWidget *parent = nullptr);
    ~E6t1();

private:
    Ui::E6t1 *ui;
};

#endif // E6T1_H
