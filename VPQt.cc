#include <Curve.h>
#include <DiscreteSpace.h>
#include <ScalarValuedFunction.h>
#include <ParametrizedFunction.h>
#include <MyTypes.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
//#include <Utilities.h>
#include <chrono>
#include <numeric>

#include <QApplication>
#include <QComboBox>
#include <QtCharts/QValueAxis>
#include <QFormLayout>
#include <QPushButton>
#include <QLabel>
#include <QDialog>
#include <QWidget>
#include <QLineEdit>
#include <QMainWindow>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <unordered_set>
#include <progressbar.hpp>
#include <RegularP1Mesh.h>
#include <RiemannLiouville.h>
#include <RiemannLiouvilleProblem.h>
#include <ParametrizedFunction.h>
#include <Msg.h>
#include <SolveRLProblem.h>

useMessages("MAIN_QT");

void plotFunction(const int &index, const BEM::CVector &vector, const std::string &fileName)
{
    QMainWindow *newWindow = new QMainWindow();
    auto functionPtr = std::unique_ptr<ScalarFunctionBase_1D>(nullptr);
    if (index == 0) {
        functionPtr = std::make_unique<PolynomialFunction_1D>(vector);
    } else if (index == 1) {
        functionPtr = std::make_unique<TrigonometricFunction_1D>(1, vector);
    } else if (index == 2) {
        functionPtr = std::make_unique<PwConstantFunction_1D>(0, 1, vector);
    }
    auto &function = *functionPtr;
    QtCharts::QLineSeries *series = new QtCharts::QLineSeries();
    double max = function(0.).real();
    double min = function(0.).real();
    for (double t = 0; t <= 1; t += 0.01) {
        double value = function(t).real();
        series->append(t, value);
        if (value > max) {
            max = value;
        } else if (value < min) {
             min = value;
        }
    }
    double padding = (max - min)*0.1;

    QtCharts::QChart *chart = new QtCharts::QChart();
    chart->addSeries(series);
    chart->createDefaultAxes();
    chart->legend()->hide();
    chart->setTitle(QString::fromStdString(fileName));
    QtCharts::QChartView *chartView = new QtCharts::QChartView(chart);
    newWindow->setCentralWidget(chartView);
    newWindow->resize(400, 300);
    newWindow->show();
}

int main(int argc, char* argv[]) {

    QApplication app(argc, argv);
    QWidget window, plotWindow;
    window.setWindowTitle("VP Qt");
    plotWindow.setWindowTitle("Plot");
    window.resize(540, 260);
    plotWindow.resize(540, 260);

    QFont font("Arial", 18);
    auto *layout = new QFormLayout(&window);
    auto *layoutPlot = new QFormLayout(&plotWindow);

    // Order Line
    auto *oRow = new QWidget(&window);
    auto *oRowLayout = new QHBoxLayout(oRow);
    oRowLayout->setContentsMargins(0, 0, 0, 0);
    oRowLayout->setAlignment(Qt::AlignVCenter);
    oRowLayout->setAlignment(Qt::AlignLeft);

    auto *orderLabel = new QLabel("Order: ");
    auto *orderEdit = new QSlider(Qt::Horizontal);
    orderEdit->setRange(51, 100);
    QLabel *orderValueLabel = new QLabel(QString::number(orderEdit->value()));
    QObject::connect(orderEdit, &QSlider::valueChanged, [&orderValueLabel](int value){
            orderValueLabel->setText(QString::number(value));});
    oRowLayout->addWidget(orderEdit);
    oRowLayout->addWidget(orderValueLabel);
    layout->addRow(orderLabel, oRow);
    
    // Difussion Line
    auto *dRow = new QWidget(&window);
    auto *dRowLayout = new QHBoxLayout(dRow);
    dRowLayout->setContentsMargins(0, 0, 0, 0);
    dRowLayout->setAlignment(Qt::AlignVCenter);
    dRowLayout->setAlignment(Qt::AlignLeft);

    auto *diffusionLabel = new QLabel("Diffusion:");
    diffusionLabel->setFont(font);
    auto *diffusionTypeBox = new QComboBox(&window);
    diffusionTypeBox->addItems({"Polynomial", "Trigonometric", "Piecewise constant"});
    auto *diffusionEdit = new QLineEdit("");
    dRowLayout->addWidget(diffusionTypeBox);
    dRowLayout->addWidget(diffusionEdit);
    auto *diffusionPlotButton = new QPushButton("Plot", &window);
    QObject::connect(diffusionPlotButton, &QPushButton::clicked, [&]() {
        QStringList lista = diffusionEdit->text().split(" ");
        BEM::CVector miVector{};
        for(const QString &s : lista) {
            miVector.push_back(s.trimmed().toDouble());
        }
        auto index = diffusionTypeBox->currentIndex();
        plotFunction(index, miVector, "Diffusion");

    }); 
    dRowLayout->addWidget(diffusionPlotButton);
    layout->addRow(diffusionLabel, dRow);

    // Reaction Line
    auto *rRow = new QWidget(&window);
    auto *rRowLayout = new QHBoxLayout(rRow);
    rRowLayout->setContentsMargins(0, 0, 0, 0);
    rRowLayout->setAlignment(Qt::AlignVCenter);
    rRowLayout->setAlignment(Qt::AlignLeft);

    auto *reactionLabel = new QLabel("Reaction:");
    auto *reactionTypeBox = new QComboBox(&window);
    reactionTypeBox->addItems({"Polynomial", "Trigonometric", "Piecewise constant"});
    auto *reactionEdit = new QLineEdit("");
    rRowLayout->addWidget(reactionTypeBox);
    rRowLayout->addWidget(reactionEdit);
    auto *reactionPlotButton = new QPushButton("Plot", &window);
    QObject::connect(reactionPlotButton, &QPushButton::clicked, [&]() {
        QStringList lista = reactionEdit->text().split(" ");
        BEM::CVector miVector{};
        for(const QString &s : lista) {
            miVector.push_back(s.trimmed().toDouble());
        }
        auto index = reactionTypeBox->currentIndex();
        plotFunction(index, miVector, "Reaction");
    });
    rRowLayout->addWidget(reactionPlotButton);
    layout->addRow(reactionLabel, rRow);

    window.show();
    //solveRLProblem(order, typeD, typeQ, dVec, qVec, rhsVec, ms);
    
    return app.exec();
}
