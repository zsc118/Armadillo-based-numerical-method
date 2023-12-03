#ifndef E4T3_H
#define E4T3_H

#include <QWidget>

namespace Ui {
class E4t3;
}

class E4t3 : public QWidget
{
    Q_OBJECT

public:
    explicit E4t3(QWidget *parent = nullptr);
    ~E4t3();

private:
    Ui::E4t3 *ui;
};

#endif // E4T3_H
