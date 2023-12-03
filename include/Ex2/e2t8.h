#ifndef E2T8_H
#define E2T8_H

#include <QWidget>

namespace Ui {
class E2t8;
}

class E2t8 : public QWidget
{
    Q_OBJECT

public:
    explicit E2t8(QWidget *parent = nullptr);
    ~E2t8();

private:
    Ui::E2t8 *ui;
};

#endif // E2T8_H
