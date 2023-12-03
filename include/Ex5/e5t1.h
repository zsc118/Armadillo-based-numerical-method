#ifndef E5T1_H
#define E5T1_H

#include <QWidget>

namespace Ui {
class E5t1;
}

class E5t1 : public QWidget
{
    Q_OBJECT

public:
    explicit E5t1(QWidget *parent = nullptr);
    ~E5t1();

private:
    Ui::E5t1 *ui;
};

#endif // E5T1_H
