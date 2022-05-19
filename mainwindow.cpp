#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "GenerationSignal.h"
#include "paramSignal.h"
#include "kol_iri_sput.h"

GenerationSignal n;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
     ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

int kol=0;

void MainWindow::on_pushButton_clicked()
{
    n.GenerateTwoLongSignal();

    for(int i = 0; i < n.listt.size(); i++)
        {
            ui->listWidget->addItem(n.listt[i]);
        }
    kol=n.listt.size();
}

void MainWindow::find_max_min(QVector<double>& x, QVector<double>y,double& maxX, double& minX, double& maxY, double& minY)
{
    for(int i=0;i<y.size();i++)
    {
        if(maxY<y[i]) maxY=y[i];
    }
    minY=-maxY/50;
    minX=x[0];
    maxX=x[x.size()-1];
}

void MainWindow::on_pushButton_2_clicked()
{
    n.CreateVFN();

    QList<QString> OneGrup;
    QList<QString> TwoGrup;
    QList<QString> ThreeGrup;
    QList<QString> SummDel;

for(int i=0;i<n.SummDelays.size();i++)
{
    OneGrup.push_back(QString(QChar::fromLatin1(n.str1[i][0])));
    TwoGrup.push_back(QString(QChar::fromLatin1(n.str1[i][2])));
    ThreeGrup.push_back(QString(QChar::fromLatin1(n.str1[i][4])));
    SummDel.push_back(QString::number(n.SummDelays[i]));
}


        // Create model:
        TestModel *PhoneBookModel = new TestModel(this);

        // Populate model with data:
        PhoneBookModel->populateData(OneGrup,TwoGrup,ThreeGrup,SummDel);

        // Connect model to table view:
        ui->tableView->setModel(PhoneBookModel);

        // Make table header visible and display table:
        ui->tableView->horizontalHeader()->setVisible(true);
        ui->tableView->show();

        for(int i=kol;i<n.listt.size();i++)
        {
            ui->listWidget->addItem(n.listt[i]);
        }
}


void MainWindow::on_listWidget_itemClicked(QListWidgetItem *item)
{
    //номер в списке
    ind=ui->listWidget->row(item);
    y.resize(n.BigMassOtrisovka[ind].size());
        x.resize(n.BigMassOtshetX[ind].size());
        for(int i=0;i<n.BigMassOtrisovka[ind].size();i++)
        {
            y[i]=n.BigMassOtrisovka[ind][i];
            x[i]=n.BigMassOtshetX[ind][i];
        }
        double maxX=0, minX=0, maxY=0, minY=0;

        find_max_min(x, y,maxX, minX,  maxY,  minY);

        // создаем график и добавляем данные:
        ui->widget->addGraph();
        ui->widget->graph(0)->setData(x, y);
        // задаем имена осей координат
        ui->widget->xAxis->setLabel("freq");
        ui->widget->yAxis->setLabel("ampl");
        // задаем размеры осей
        ui->widget->xAxis->setRange(minX, maxX);
        ui->widget->yAxis->setRange(minY, maxY);
        ui->widget->replot();
}


TestModel::TestModel(QObject *parent) : QAbstractTableModel(parent)
{
}

// Create a method to populate the model with data:
void TestModel::populateData(const QList<QString> &OneGrup,const QList<QString> &TwoGrup, const QList<QString> &ThreeGrup,const QList<QString>& SummDel)
{
    tm_one_grup.clear();
    tm_one_grup = OneGrup;
    tm_two_grup.clear();
    tm_two_grup = TwoGrup;
    tm_three_grup.clear();
    tm_three_grup = ThreeGrup;
    tm_summ_del.clear();
    tm_summ_del = SummDel;
    return;
}

int TestModel::rowCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent);
    return tm_one_grup.length();
}

int TestModel::columnCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent);
    return 4;
}

QVariant TestModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid() || role != Qt::DisplayRole) {
        return QVariant();
    }
    if (index.column() == 0)
    {
        return tm_one_grup[index.row()];
    } else if (index.column() == 1)
    {
        return tm_two_grup[index.row()];
    } else if (index.column() == 2)
    {
        return tm_three_grup[index.row()];
    }
    else if (index.column() == 3)
        {
            return tm_summ_del[index.row()];
        }
    return QVariant();
}

QVariant TestModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole && orientation == Qt::Horizontal)
    {
        if (section == 0) return QString("1 спутник");
        else if (section == 1) return QString("2 спутник");
        else if (section == 2) return QString("3 спутник");
        else if (section == 3) return QString("Итог критерия");
    }
    return QVariant();
}

void MainWindow::on_action_2_triggered()
{
    Dialog *frm = new Dialog();
    frm->show();
}


void MainWindow::on_action_5_triggered()
{
    Kol_IRI_Sput *frm = new Kol_IRI_Sput();
    frm->show();
}

