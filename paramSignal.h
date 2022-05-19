#ifndef PARAMSIGNAL_H
#define PARAMSIGNAL_H

#include <QDialog>
#include "mainwindow.h"

namespace Ui {
class Dialog;
}

class TableParam : public QAbstractTableModel
{
    Q_OBJECT

public:
    TableParam(QObject *parent = 0);

    void populateData(const QList<QString> &Name,const QList<QString> &Znach);

    int rowCount(const QModelIndex &parent = QModelIndex()) const Q_DECL_OVERRIDE;
    int columnCount(const QModelIndex &parent = QModelIndex()) const Q_DECL_OVERRIDE;

    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const Q_DECL_OVERRIDE;
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const Q_DECL_OVERRIDE;

private:
    QList<QString> tm_name;
    QList<QString> tm_znach;
};

class Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit Dialog(QWidget *parent = nullptr);
    ~Dialog();




private:
    Ui::Dialog *ui;
};

#endif // PARAMSIGNAL_H

