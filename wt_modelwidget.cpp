/*
 * 文件名: wt_modelwidget.cpp
 * 文件作用: 压裂水平井复合页岩油模型界面类实现
 * 功能描述:
 * 1. 界面初始化：根据模型类型动态调整参数输入框的显示/隐藏。
 * 2. 求解器集成：根据 ID 范围分发计算任务。
 * 3. 业务逻辑：处理计算请求、参数读取、结果绘图与导出。
 * 修改记录:
 * 1. [修复] 移除动态创建代码，改为直接使用 UI 文件中的 alphaEdit 和 cPhiEdit 控件。
 * 2. [修复] 修正了参数可见性逻辑，确保切换模型时控件状态正确。
 */

#include "wt_modelwidget.h"
#include "ui_wt_modelwidget.h"
#include "modelmanager.h"
#include "modelparameter.h"

#include <QDebug>
#include <QMessageBox>
#include <QFileDialog>
#include <QTextStream>
#include <QDateTime>
#include <QCoreApplication>
#include <QSplitter>
#include <QLabel>
#include <QLineEdit>
#include <QGridLayout>
#include <cmath>

WT_ModelWidget::WT_ModelWidget(ModelType type, QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::WT_ModelWidget)
    , m_type(type)
    , m_solver1(nullptr)
    , m_solver2(nullptr)
    , m_highPrecision(true)
{
    ui->setupUi(this);

    // 实例化求解器
    // 0-35: Group 1 (Solver1)
    // 36-71: Group 2 (Solver2)
    if (m_type >= 0 && m_type <= 35) {
        m_solver1 = new ModelSolver01_06((ModelSolver01_06::ModelType)m_type);
    }
    else if (m_type >= 36 && m_type <= 71) {
        m_solver2 = new ModelSolver19_36((ModelSolver19_36::ModelType)(m_type - 36));
    }

    m_colorList = { Qt::red, Qt::blue, QColor(0,180,0), Qt::magenta, QColor(255,140,0), Qt::cyan };

    QList<int> sizes;
    sizes << 240 << 960;
    ui->splitter->setSizes(sizes);
    ui->splitter->setCollapsible(0, false);

    ui->btnSelectModel->setText(getModelName());

    initUi();
    initChart();
    setupConnections();
    onResetParameters();
}

WT_ModelWidget::~WT_ModelWidget()
{
    if(m_solver1) delete m_solver1;
    if(m_solver2) delete m_solver2;
    delete ui;
}

QString WT_ModelWidget::getModelName() const {
    if (m_solver1) return ModelSolver01_06::getModelName((ModelSolver01_06::ModelType)m_type, false);
    if (m_solver2) return ModelSolver19_36::getModelName((ModelSolver19_36::ModelType)(m_type - 36), false);
    return "未知模型";
}

WT_ModelWidget::ModelCurveData WT_ModelWidget::calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime)
{
    if (m_solver1) return m_solver1->calculateTheoreticalCurve(params, providedTime);
    if (m_solver2) return m_solver2->calculateTheoreticalCurve(params, providedTime);
    return ModelCurveData();
}

void WT_ModelWidget::setHighPrecision(bool high)
{
    m_highPrecision = high;
    if (m_solver1) m_solver1->setHighPrecision(high);
    if (m_solver2) m_solver2->setHighPrecision(high);
}

void WT_ModelWidget::initUi() {
    // 1. 基础标签更新
    if(ui->label_km) ui->label_km->setText("流度比 M12");
    if(ui->label_rmD) ui->label_rmD->setText("复合半径 rm (m)");
    if(ui->label_reD) ui->label_reD->setText("外区半径 re (m)");
    if(ui->label_cD) ui->label_cD->setText("井筒储集 C (m³/MPa)");

    // --- 2. 可见性控制逻辑 ---

    // A. 边界条件 (每12个一组，0-3为无限大)
    int groupIdx = m_type % 12;
    bool isInfinite = (groupIdx < 4);
    ui->label_reD->setVisible(!isInfinite);
    ui->reDEdit->setVisible(!isInfinite);

    // B. 井储与表皮逻辑
    // 0: Const, 1: Line, 2: Fair, 3: Hegeman
    int storageType = m_type % 4;

    // 线源解(1)不显示 C 和 S，其他(0, 2, 3)都需要显示
    bool hasBasicStorage = (storageType != 1);
    ui->label_cD->setVisible(hasBasicStorage);
    ui->cDEdit->setVisible(hasBasicStorage);
    ui->label_s->setVisible(hasBasicStorage);
    ui->sEdit->setVisible(hasBasicStorage);

    // [新增] 变井储参数可见性
    // 仅 Fair(2) 和 Hegeman(3) 需要 alpha 和 C_phi
    bool hasVarStorage = (storageType == 2 || storageType == 3);

    ui->label_alpha->setVisible(hasVarStorage);
    ui->alphaEdit->setVisible(hasVarStorage);
    ui->label_cphi->setVisible(hasVarStorage);
    ui->cPhiEdit->setVisible(hasVarStorage);

    // C. 介质参数可见性
    bool hasInnerParams = false;
    bool hasOuterParams = false;

    if (m_type <= 35) { // Solver 1
        if (m_type <= 23) hasInnerParams = true; // 0-23 Inner Dual
        if (m_type <= 11) hasOuterParams = true; // 0-11 Outer Dual
    } else { // Solver 2
        hasInnerParams = true; // 页岩型内区始终需要
        int sub = m_type - 36;
        if (sub <= 11 || sub >= 24) hasOuterParams = true; // Shale or Dual outer
    }

    ui->label_omga1->setVisible(hasInnerParams);
    ui->omga1Edit->setVisible(hasInnerParams);
    ui->label_remda1->setVisible(hasInnerParams);
    ui->remda1Edit->setVisible(hasInnerParams);

    ui->label_omga2->setVisible(hasOuterParams);
    ui->omga2Edit->setVisible(hasOuterParams);

    ui->label_remda2->setVisible(hasOuterParams);
    ui->remda2Edit->setVisible(hasOuterParams);
}

void WT_ModelWidget::initChart() {
    MouseZoom* plot = ui->chartWidget->getPlot();
    plot->setBackground(Qt::white);
    plot->axisRect()->setBackground(Qt::white);

    QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
    plot->xAxis->setScaleType(QCPAxis::stLogarithmic); plot->xAxis->setTicker(logTicker);
    plot->yAxis->setScaleType(QCPAxis::stLogarithmic); plot->yAxis->setTicker(logTicker);
    plot->xAxis->setNumberFormat("eb"); plot->xAxis->setNumberPrecision(0);
    plot->yAxis->setNumberFormat("eb"); plot->yAxis->setNumberPrecision(0);

    QFont labelFont("Microsoft YaHei", 10, QFont::Bold);
    QFont tickFont("Microsoft YaHei", 9);
    plot->xAxis->setLabel("时间 Time (h)");
    plot->yAxis->setLabel("压力 & 导数 Pressure & Derivative (MPa)");
    plot->xAxis->setLabelFont(labelFont); plot->yAxis->setLabelFont(labelFont);
    plot->xAxis->setTickLabelFont(tickFont); plot->yAxis->setTickLabelFont(tickFont);

    plot->xAxis2->setVisible(true); plot->yAxis2->setVisible(true);
    plot->xAxis2->setTickLabels(false); plot->yAxis2->setTickLabels(false);
    connect(plot->xAxis, SIGNAL(rangeChanged(QCPRange)), plot->xAxis2, SLOT(setRange(QCPRange)));
    connect(plot->yAxis, SIGNAL(rangeChanged(QCPRange)), plot->yAxis2, SLOT(setRange(QCPRange)));
    plot->xAxis2->setScaleType(QCPAxis::stLogarithmic); plot->yAxis2->setScaleType(QCPAxis::stLogarithmic);
    plot->xAxis2->setTicker(logTicker); plot->yAxis2->setTicker(logTicker);

    plot->xAxis->grid()->setVisible(true); plot->yAxis->grid()->setVisible(true);
    plot->xAxis->grid()->setSubGridVisible(true); plot->yAxis->grid()->setSubGridVisible(true);
    plot->xAxis->grid()->setPen(QPen(QColor(220, 220, 220), 1, Qt::SolidLine));
    plot->yAxis->grid()->setPen(QPen(QColor(220, 220, 220), 1, Qt::SolidLine));
    plot->xAxis->grid()->setSubGridPen(QPen(QColor(240, 240, 240), 1, Qt::DotLine));
    plot->yAxis->grid()->setSubGridPen(QPen(QColor(240, 240, 240), 1, Qt::DotLine));

    plot->xAxis->setRange(1e-3, 1e3); plot->yAxis->setRange(1e-3, 1e2);

    plot->legend->setVisible(true);
    plot->legend->setFont(QFont("Microsoft YaHei", 9));
    plot->legend->setBrush(QBrush(QColor(255, 255, 255, 200)));

    ui->chartWidget->setTitle("复合页岩油储层试井曲线");
}

void WT_ModelWidget::setupConnections() {
    connect(ui->calculateButton, &QPushButton::clicked, this, &WT_ModelWidget::onCalculateClicked);
    connect(ui->resetButton, &QPushButton::clicked, this, &WT_ModelWidget::onResetParameters);
    connect(ui->chartWidget, &ChartWidget::exportDataTriggered, this, &WT_ModelWidget::onExportData);
    connect(ui->btnExportDataTab, &QPushButton::clicked, this, &WT_ModelWidget::onExportData);
    connect(ui->checkShowPoints, &QCheckBox::toggled, this, &WT_ModelWidget::onShowPointsToggled);
    connect(ui->btnSelectModel, &QPushButton::clicked, this, &WT_ModelWidget::requestModelSelection);
}

QVector<double> WT_ModelWidget::parseInput(const QString& text) {
    QVector<double> values;
    QString cleanText = text;
    cleanText.replace("，", ",");
    QStringList parts = cleanText.split(",", Qt::SkipEmptyParts);
    for(const QString& part : parts) {
        bool ok;
        double v = part.trimmed().toDouble(&ok);
        if(ok) values.append(v);
    }
    if(values.isEmpty()) values.append(0.0);
    return values;
}

void WT_ModelWidget::setInputText(QLineEdit* edit, double value) {
    if(!edit) return;
    edit->setText(QString::number(value, 'g', 8));
}

void WT_ModelWidget::onResetParameters() {
    ModelManager mgr;
    QMap<QString, double> defaults = mgr.getDefaultParameters(static_cast<ModelManager::ModelType>(m_type));

    // 基础参数
    setInputText(ui->phiEdit, defaults["phi"]);
    setInputText(ui->hEdit, defaults["h"]);
    setInputText(ui->rwEdit, defaults["rw"]);
    setInputText(ui->muEdit, defaults["mu"]);
    setInputText(ui->BEdit, defaults["B"]);
    setInputText(ui->CtEdit, defaults["Ct"]);
    setInputText(ui->qEdit, defaults["q"]);
    setInputText(ui->tEdit, 1000.0);
    setInputText(ui->pointsEdit, 100);

    // 模型参数
    setInputText(ui->kfEdit, defaults["kf"]);
    setInputText(ui->kmEdit, defaults["M12"]);
    setInputText(ui->LEdit, defaults["L"]);
    setInputText(ui->LfEdit, defaults["Lf"]);
    setInputText(ui->nfEdit, defaults["nf"]);
    setInputText(ui->rmDEdit, defaults["rm"]);

    setInputText(ui->omga1Edit, defaults.value("omega1", 0.0));
    setInputText(ui->omga2Edit, defaults.value("omega2", 0.0));
    setInputText(ui->remda1Edit, defaults.value("lambda1", 0.0));
    setInputText(ui->gamaDEdit, defaults.value("gamaD", 0.02));

    setInputText(ui->remda2Edit, defaults.value("lambda2", 0.0));
    setInputText(ui->eta12Edit, defaults.value("eta12", 0.2));

    // [新增] 变井储参数回显
    setInputText(ui->alphaEdit, defaults.value("alpha", 0.1));
    setInputText(ui->cPhiEdit, defaults.value("C_phi", 1e-4));

    if (ui->reDEdit->isVisible()) {
        setInputText(ui->reDEdit, defaults.value("re", 20000.0));
    }

    if (ui->cDEdit->isVisible()) {
        setInputText(ui->cDEdit, 0.1);
        setInputText(ui->sEdit, defaults.value("S", 0.0));
    }
}

void WT_ModelWidget::onDependentParamsChanged() {}

void WT_ModelWidget::onShowPointsToggled(bool checked) {
    MouseZoom* plot = ui->chartWidget->getPlot();
    for(int i = 0; i < plot->graphCount(); ++i) {
        if (checked) plot->graph(i)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 5));
        else plot->graph(i)->setScatterStyle(QCPScatterStyle::ssNone);
    }
    plot->replot();
}

void WT_ModelWidget::onCalculateClicked() {
    ui->calculateButton->setEnabled(false);
    ui->calculateButton->setText("计算中...");
    QCoreApplication::processEvents();
    runCalculation();
    ui->calculateButton->setEnabled(true);
    ui->calculateButton->setText("开始计算");
}

void WT_ModelWidget::runCalculation() {
    MouseZoom* plot = ui->chartWidget->getPlot();
    plot->clearGraphs();

    QMap<QString, QVector<double>> rawParams;

    // 读取界面参数
    rawParams["phi"] = parseInput(ui->phiEdit->text());
    rawParams["h"] = parseInput(ui->hEdit->text());
    rawParams["rw"] = parseInput(ui->rwEdit->text());
    rawParams["mu"] = parseInput(ui->muEdit->text());
    rawParams["B"] = parseInput(ui->BEdit->text());
    rawParams["Ct"] = parseInput(ui->CtEdit->text());
    rawParams["q"] = parseInput(ui->qEdit->text());
    rawParams["t"] = parseInput(ui->tEdit->text());

    rawParams["kf"] = parseInput(ui->kfEdit->text());
    rawParams["M12"] = parseInput(ui->kmEdit->text());
    rawParams["L"] = parseInput(ui->LEdit->text());
    rawParams["Lf"] = parseInput(ui->LfEdit->text());
    rawParams["nf"] = parseInput(ui->nfEdit->text());
    rawParams["rm"] = parseInput(ui->rmDEdit->text());
    rawParams["omega1"] = parseInput(ui->omga1Edit->text());
    rawParams["omega2"] = parseInput(ui->omga2Edit->text());
    rawParams["lambda1"] = parseInput(ui->remda1Edit->text());
    rawParams["gamaD"] = parseInput(ui->gamaDEdit->text());

    rawParams["lambda2"] = parseInput(ui->remda2Edit->text());
    rawParams["eta12"] = parseInput(ui->eta12Edit->text());

    // [新增] 读取变井储参数
    if (ui->alphaEdit->isVisible()) rawParams["alpha"] = parseInput(ui->alphaEdit->text());
    else rawParams["alpha"] = {0.1};

    if (ui->cPhiEdit->isVisible()) rawParams["C_phi"] = parseInput(ui->cPhiEdit->text());
    else rawParams["C_phi"] = {1e-4};

    if (ui->reDEdit->isVisible()) rawParams["re"] = parseInput(ui->reDEdit->text());
    else rawParams["re"] = {20000.0};

    // 井储系数 C -> cD 转换
    if (ui->cDEdit->isVisible()) {
        QVector<double> C_vals = parseInput(ui->cDEdit->text());
        QVector<double> cD_vals;

        double phi = rawParams["phi"].isEmpty() ? 0.05 : rawParams["phi"].first();
        double Ct = rawParams["Ct"].isEmpty() ? 5e-4 : rawParams["Ct"].first();
        double h = rawParams["h"].isEmpty() ? 20.0 : rawParams["h"].first();
        double L = rawParams["L"].isEmpty() ? 1000.0 : rawParams["L"].first();

        double denom = phi * h * Ct * L * L;
        double factor = (denom > 1e-20) ? 0.159 / denom : 0.0;

        for(double valC : C_vals) cD_vals.append(valC * factor);

        rawParams["cD"] = cD_vals;
        rawParams["S"] = parseInput(ui->sEdit->text());
    } else {
        rawParams["cD"] = {0.0};
        rawParams["S"] = {0.0};
    }

    // 敏感性分析逻辑
    QString sensitivityKey = "";
    QVector<double> sensitivityValues;
    for(auto it = rawParams.begin(); it != rawParams.end(); ++it) {
        if(it.key() == "t") continue;
        if(it.value().size() > 1) {
            sensitivityKey = it.key();
            sensitivityValues = it.value();
            break;
        }
    }
    bool isSensitivity = !sensitivityKey.isEmpty();

    QMap<QString, double> baseParams;
    for(auto it = rawParams.begin(); it != rawParams.end(); ++it) {
        baseParams[it.key()] = it.value().isEmpty() ? 0.0 : it.value().first();
    }
    baseParams["N"] = m_highPrecision ? 10.0 : 4.0;
    if(baseParams["L"] > 1e-9) baseParams["LfD"] = baseParams["Lf"] / baseParams["L"];

    int nPoints = ui->pointsEdit->text().toInt();
    if(nPoints < 5) nPoints = 5;
    double maxTime = baseParams.value("t", 1000.0);
    QVector<double> t = ModelManager::generateLogTimeSteps(nPoints, -3.0, log10(maxTime));

    int iterations = isSensitivity ? sensitivityValues.size() : 1;
    iterations = qMin(iterations, (int)m_colorList.size());

    QString resultTextHeader = QString("计算完成 (%1)\n").arg(getModelName());
    if(isSensitivity) resultTextHeader += QString("敏感性参数: %1\n").arg(sensitivityKey);

    for(int i = 0; i < iterations; ++i) {
        QMap<QString, double> currentParams = baseParams;
        double val = 0;
        if (isSensitivity) {
            val = sensitivityValues[i];
            currentParams[sensitivityKey] = val;
            if (sensitivityKey == "L" || sensitivityKey == "Lf") {
                if(currentParams["L"] > 1e-9) currentParams["LfD"] = currentParams["Lf"] / currentParams["L"];
            }
        }

        ModelCurveData res = calculateTheoreticalCurve(currentParams, t);
        res_tD = std::get<0>(res);
        res_pD = std::get<1>(res);
        res_dpD = std::get<2>(res);

        QColor curveColor = isSensitivity ? m_colorList[i] : Qt::red;
        QString legendName = isSensitivity ? QString("%1 = %2").arg(sensitivityKey).arg(val) : "理论曲线";
        plotCurve(res, legendName, curveColor, isSensitivity);
    }

    QString resultText = resultTextHeader;
    resultText += "t(h)\t\tDp(MPa)\t\tdDp(MPa)\n";
    for(int i=0; i<res_pD.size(); ++i) {
        resultText += QString("%1\t%2\t%3\n").arg(res_tD[i],0,'e',4).arg(res_pD[i],0,'e',4).arg(res_dpD[i],0,'e',4);
    }
    ui->resultTextEdit->setText(resultText);

    ui->chartWidget->getPlot()->rescaleAxes();
    plot->replot();
    onShowPointsToggled(ui->checkShowPoints->isChecked());
    emit calculationCompleted(getModelName(), baseParams);
}

void WT_ModelWidget::plotCurve(const ModelCurveData& data, const QString& name, QColor color, bool isSensitivity) {
    MouseZoom* plot = ui->chartWidget->getPlot();
    const QVector<double>& t = std::get<0>(data);
    const QVector<double>& p = std::get<1>(data);
    const QVector<double>& d = std::get<2>(data);

    QCPGraph* graphP = plot->addGraph();
    graphP->setData(t, p);
    graphP->setPen(QPen(color, 2, Qt::SolidLine));

    QCPGraph* graphD = plot->addGraph();
    graphD->setData(t, d);

    if (isSensitivity) {
        graphD->setPen(QPen(color, 2, Qt::DashLine));
        graphP->setName(name);
        graphD->removeFromLegend();
    } else {
        graphP->setPen(QPen(Qt::red, 2));
        graphP->setName("压力");
        graphD->setPen(QPen(Qt::blue, 2));
        graphD->setName("压力导数");
    }
}

void WT_ModelWidget::onExportData() {
    if (res_tD.isEmpty()) return;
    QString defaultDir = ModelParameter::instance()->getProjectPath();
    if(defaultDir.isEmpty()) defaultDir = ".";
    QString path = QFileDialog::getSaveFileName(this, "导出CSV数据", defaultDir + "/CalculatedData.csv", "CSV Files (*.csv)");
    if (path.isEmpty()) return;
    QFile f(path);
    if (f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&f);
        out << "t,Dp,dDp\n";
        for (int i = 0; i < res_tD.size(); ++i) {
            double dp = (i < res_dpD.size()) ? res_dpD[i] : 0.0;
            out << res_tD[i] << "," << res_pD[i] << "," << dp << "\n";
        }
        f.close();
        QMessageBox::information(this, "导出成功", "数据文件已保存");
    }
}
