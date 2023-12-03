#ifndef E1T2_H
#define E1T2_H

#include <QWidget>
#include <string>
#include <vector>
namespace Ui {
class e1t2;
}

class e1t2 : public QWidget
{
    Q_OBJECT

public:
    explicit e1t2(QWidget *parent = nullptr);
    ~e1t2();

protected:
    std::vector<std::string> func;// 迭代函数
    std::vector<std::string> x;// 未知数
    std::vector<double> current_vec;// 当前迭代向量值
    Ui::e1t2 *ui;
};

#endif // E1T2_H
