#ifndef E5T2_H
#define E5T2_H

#include <QWidget>

namespace Ui {
class E5t2;
}

class E5t2 : public QWidget
{
    Q_OBJECT

public:
    explicit E5t2(QWidget *parent = nullptr);
    ~E5t2();

private:
    Ui::E5t2 *ui;
};

#endif // E5T2_H
