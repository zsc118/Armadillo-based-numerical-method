#ifndef E4T4_H
#define E4T4_H

#include <QWidget>

namespace Ui {
class E4t4;
}

class E4t4 : public QWidget
{
    Q_OBJECT

public:
    explicit E4t4(QWidget *parent = nullptr);
    ~E4t4();

private:
    Ui::E4t4 *ui;
};

#endif // E4T4_H
