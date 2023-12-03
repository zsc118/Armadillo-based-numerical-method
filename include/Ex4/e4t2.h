#ifndef E4T2_H
#define E4T2_H

#include <QWidget>

namespace Ui {
class E4t2;
}

class E4t2 : public QWidget
{
    Q_OBJECT

public:
    explicit E4t2(QWidget *parent = nullptr);
    ~E4t2();

private:
    Ui::E4t2 *ui;
};

#endif // E4T2_H
