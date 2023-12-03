#ifndef E5T3_H
#define E5T3_H

#include <QWidget>

namespace Ui {
class E5t3;
}

class E5t3 : public QWidget
{
    Q_OBJECT

public:
    explicit E5t3(QWidget *parent = nullptr);
    ~E5t3();

private:
    Ui::E5t3 *ui;
};

#endif // E5T3_H
