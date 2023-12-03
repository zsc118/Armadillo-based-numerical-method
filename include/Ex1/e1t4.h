#ifndef E1T4_H
#define E1T4_H

#include <QWidget>

namespace Ui {
class E1t4;
}

class E1t4 : public QWidget
{
    Q_OBJECT

public:
    explicit E1t4(QWidget *parent = nullptr);
    ~E1t4();

private:
    Ui::E1t4 *ui;
};

#endif // E1T4_H
