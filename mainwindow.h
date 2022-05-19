#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"
#include <math.h>
#include <vector>
using namespace std;

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class TestModel : public QAbstractTableModel
{
    Q_OBJECT

public:
    TestModel(QObject *parent = 0);

    void populateData(const QList<QString> &OneGrup,const QList<QString> &TwoGrup, const QList<QString> &ThreeGrup, const QList<QString> &SummDel);

    int rowCount(const QModelIndex &parent = QModelIndex()) const Q_DECL_OVERRIDE;
    int columnCount(const QModelIndex &parent = QModelIndex()) const Q_DECL_OVERRIDE;

    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const Q_DECL_OVERRIDE;
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const Q_DECL_OVERRIDE;

private:
    QList<QString> tm_one_grup;
    QList<QString> tm_two_grup;
    QList<QString> tm_three_grup;
    QList<QString> tm_summ_del;
};



class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
private slots:
     void on_pushButton_clicked();
     void on_pushButton_2_clicked();
     void on_listWidget_itemClicked(QListWidgetItem *item);
     void find_max_min(QVector<double>& x, QVector<double>y,double& maxX, double& minX, double& maxY, double& minY);

     void on_action_2_triggered();

     void on_action_5_triggered();

public:
     QVector <double>x,y;
     int ind;
     double bufmaxY;
};
#endif // MAINWINDOW_H
