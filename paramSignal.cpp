#include "paramSignal.h"
#include "ui_paramSignal.h"
#include "GenerationSignal.h"

GenerationSignal s;

Dialog::Dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog)
{
    ui->setupUi(this);
    QList<QString> Name;
    QList<QString> Znach;



    Name.push_back(QString::fromStdString("Частота дискретизации"));
    Name.push_back(QString::fromStdString("Метка времени начала"));
    Name.push_back(QString::fromStdString("Продолжительность"));
    Name.push_back(QString::fromStdString("Начальная фаза"));
    Name.push_back(QString::fromStdString("Нечто"));
    Name.push_back(QString::fromStdString("Скорость передачи данных"));
    Name.push_back(QString::fromStdString("Дополнительный параметр"));

    Znach.push_back(QString::number(s.samplingFrequency));
    Znach.push_back(QString::number(s.startTimestamp));
    Znach.push_back(QString::number(s.Duration));
    Znach.push_back(QString::number(s.startPhase));
    Znach.push_back(QString::number(s.nSamples));
    Znach.push_back(QString::number(s.Bitrate));
    Znach.push_back(QString::number(s.additionalParameter));


        // Create model:
        TableParam *Param = new TableParam(this);

        // Populate model with data:
        Param->populateData(Name,Znach);

        // Connect model to table view:
        ui->tableView->setModel(Param);

        // Make table header visible and display table:
        ui->tableView->horizontalHeader()->setVisible(true);
        ui->tableView->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
        ui->tableView->show();
}

Dialog::~Dialog()
{
    delete ui;
}
TableParam::TableParam(QObject *parent) : QAbstractTableModel(parent)
{
}

// Create a method to populate the model with data:
void TableParam::populateData(const QList<QString> &Name,const QList<QString> &Znach)
{
    tm_name.clear();
    tm_name = Name;
    tm_znach.clear();
    tm_znach = Znach;
    return;
}

int TableParam::rowCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent);
    return tm_name.length();
}

int TableParam::columnCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent);
    return 2;
}

QVariant TableParam::data(const QModelIndex &index, int role) const
{
    if (!index.isValid() || role != Qt::DisplayRole) {
        return QVariant();
    }
    if (index.column() == 0)
    {
        return tm_name[index.row()];
    } else if (index.column() == 1)
    {
        return tm_znach[index.row()];
    }
    return QVariant();
}

QVariant TableParam::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole && orientation == Qt::Horizontal)
    {
        if (section == 0) return QString("Название");
        else if (section == 1) return QString("Значение");
    }
    return QVariant();
}



