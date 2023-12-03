#ifndef E2T4_H
#define E2T4_H

#include <QWidget>

namespace Ui {
class E2t4;
}

class E2t4 : public QWidget
{
    Q_OBJECT

public:
    explicit E2t4(QWidget *parent = nullptr);
    ~E2t4();

private:
    Ui::E2t4 *ui;
};

#endif // E2T4_H
