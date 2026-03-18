/*
 * 文件名: paramselectdialog.cpp
 * 文件作用: 参数选择配置对话框的具体实现文件。
 * 功能描述:
 * 1. 动态生成参数配置表格：实现了“拟合变量”第一列、“显示”最后一列的全新列顺序。
 * 2. 完美的UI细节：使用“透明控件+底层背景色”方案，彻底解决了高亮色遮挡单元格边框线的问题。
 * 3. 优化高亮体验：高亮采用淡黄色，且仅高亮中间数据区域，首尾复选框列保持默认交替色。
 * 4. 列表文字内容上下、左右完全居中，行高增加以提升阅读体验。
 * 5. 实现了数据的双向收集与自动上下限运算。
 */

#include "paramselectdialog.h"
#include "ui_paramselectdialog.h"
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QDebug>
#include <QMessageBox>

// 自定义 SpinBox：重写 textFromValue 去除末尾多余的 0
class SmartDoubleSpinBox : public QDoubleSpinBox {
public:
    explicit SmartDoubleSpinBox(QWidget* parent = nullptr) : QDoubleSpinBox(parent) {}
    QString textFromValue(double value) const override {
        return QString::number(value, 'g', decimals());
    }
};

// 构造函数
ParamSelectDialog::ParamSelectDialog(const QList<FitParameter> &params,
                                     ModelManager::ModelType modelType,
                                     double fitTime,
                                     QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ParamSelectDialog),
    m_params(params),
    m_modelType(modelType)
{
    ui->setupUi(this);
    this->setWindowTitle("拟合参数配置");

    ui->spinTimeMax->setValue(fitTime);

    connect(ui->btnOk, &QPushButton::clicked, this, &ParamSelectDialog::onConfirm);
    connect(ui->btnCancel, &QPushButton::clicked, this, &ParamSelectDialog::onCancel);
    connect(ui->btnResetDefaults, &QPushButton::clicked, this, &ParamSelectDialog::onResetParams);
    connect(ui->btnAutoLimits, &QPushButton::clicked, this, &ParamSelectDialog::onAutoLimits);

    ui->btnCancel->setAutoDefault(false);

    // 初始化表格
    initTable();
}

// 析构函数
ParamSelectDialog::~ParamSelectDialog()
{
    delete ui;
}

// 拦截滚轮事件
bool ParamSelectDialog::eventFilter(QObject *obj, QEvent *event)
{
    if (event->type() == QEvent::Wheel) {
        if (qobject_cast<QAbstractSpinBox*>(obj)) return true;
    }
    return QDialog::eventFilter(obj, event);
}

// 槽函数：重置默认参数
void ParamSelectDialog::onResetParams()
{
    if (QMessageBox::question(this, "确认", "确定要重置为该模型的默认参数吗？") != QMessageBox::Yes) return;
    m_params = FittingParameterChart::generateDefaultParams(m_modelType);
    FittingParameterChart::adjustLimits(m_params);
    initTable();
}

// 槽函数：自动更新限值
void ParamSelectDialog::onAutoLimits()
{
    collectData();
    FittingParameterChart::adjustLimits(m_params);
    initTable();
    QMessageBox::information(this, "提示", "参数上下限及滚轮步长已根据当前值更新。");
}

// 核心渲染函数：更新行高亮状态（首尾不亮，纯粹底层高亮以显示边框）
void ParamSelectDialog::updateRowAppearance(int row, bool isFit)
{
    QString yellowHex = "#FFFFCC"; // 淡黄色

    // 遍历该行的所有列
    for(int col = 0; col < ui->tableWidget->columnCount(); ++col) {
        QTableWidgetItem* item = ui->tableWidget->item(row, col);
        if(item) {
            // [优化] 只高亮中间区域，第0列(拟合变量)和第7列(显示)保持原样交替色
            bool shouldHighlight = isFit && (col != 0 && col != 7);
            if(shouldHighlight) {
                item->setBackground(QColor(yellowHex));
            } else {
                item->setData(Qt::BackgroundRole, QVariant()); // 取消自定义背景
            }
        }
        // 注意：内部所有嵌在单元格内的控件（SpinBox, Widget）在 initTable 中
        // 已经被强制设定为 transparent（透明），因此底层的颜色会透出来，并且绝对不会遮挡外层网格线。
    }
}

// 内部函数：初始化表格视图
void ParamSelectDialog::initTable()
{
    ui->tableWidget->clear();
    QStringList headers;
    // [重排] 更新列表头的顺序：拟合变量(首) -> 显示(尾)
    headers << "拟合变量" << "当前数值" << "单位" << "参数名称" << "下限" << "上限" << "滚轮步长" << "显示";
    ui->tableWidget->setColumnCount(headers.size());
    ui->tableWidget->setHorizontalHeaderLabels(headers);
    ui->tableWidget->setRowCount(m_params.size());

    // 整体表格美化：交替行色、表头样式、默认行高加大
    ui->tableWidget->setAlternatingRowColors(true);
    ui->tableWidget->horizontalHeader()->setStyleSheet("QHeaderView::section { background-color: #f0f0f0; border: 1px solid #dcdcdc; padding: 4px; font-weight: bold; }");
    ui->tableWidget->verticalHeader()->setVisible(false);
    ui->tableWidget->verticalHeader()->setDefaultSectionSize(40); // 增加行高

    // 设置表格基本网格颜色及内部 Item 内边距
    ui->tableWidget->setStyleSheet("QTableWidget { border: 1px solid #dcdcdc; gridline-color: #e0e0e0; selection-background-color: transparent; } "
                                   "QTableWidget::item { padding: 4px; }");

    QString checkBoxStyle =
        "QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #cccccc; border-radius: 3px; background-color: white; }"
        "QCheckBox::indicator:checked { background-color: #0078d7; border-color: #0078d7; }"
        "QCheckBox::indicator:hover { border-color: #0078d7; }";

    // 一次性定义透明控件样式，防止遮挡边框
    QString transparentSpinStyle = "QDoubleSpinBox { background-color: transparent; border: none; margin: 0px; }";
    QString transparentWidgetStyle = "QWidget { background-color: transparent; }";

    for(int i = 0; i < m_params.size(); ++i) {
        const FitParameter& p = m_params[i];

        // 预分配底层的 Item，用于承载高亮颜色
        for (int col = 0; col < 8; ++col) {
            if (!ui->tableWidget->item(i, col)) {
                ui->tableWidget->setItem(i, col, new QTableWidgetItem());
            }
        }

        // --- 第 0 列: 拟合变量 (复选框) ---
        QWidget* pWidgetFit = new QWidget();
        pWidgetFit->setStyleSheet(transparentWidgetStyle);
        QHBoxLayout* pLayoutFit = new QHBoxLayout(pWidgetFit);
        QCheckBox* chkFit = new QCheckBox();
        chkFit->setChecked(p.isFit);
        chkFit->setStyleSheet(checkBoxStyle);
        pLayoutFit->addWidget(chkFit);
        pLayoutFit->setAlignment(Qt::AlignCenter);
        pLayoutFit->setContentsMargins(0,0,0,0);
        if (p.name == "LfD") { chkFit->setEnabled(false); chkFit->setChecked(false); }
        ui->tableWidget->setCellWidget(i, 0, pWidgetFit);

        // --- 第 1 列: 当前数值 (调节框) ---
        SmartDoubleSpinBox* spinVal = new SmartDoubleSpinBox();
        spinVal->setStyleSheet(transparentSpinStyle); // 完全透明不遮挡
        spinVal->setRange(-9e9, 9e9);
        spinVal->setDecimals(10);
        spinVal->setValue(p.value);
        spinVal->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 1, spinVal);

        // --- 第 2 列: 单位 (文本) ---
        QString dummy, dummy2, dummy3, unitStr;
        FittingParameterChart::getParamDisplayInfo(p.name, dummy, dummy2, dummy3, unitStr);
        if(unitStr == "无因次" || unitStr == "小数") unitStr = "-";
        QTableWidgetItem* unitItem = ui->tableWidget->item(i, 2);
        unitItem->setText(unitStr);
        unitItem->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter); // 完全居中
        unitItem->setFlags(unitItem->flags() & ~Qt::ItemIsEditable);

        // --- 第 3 列: 参数名称 (文本) ---
        QString displayNameFull = QString("%1 (%2)").arg(p.displayName).arg(p.name);
        QTableWidgetItem* nameItem = ui->tableWidget->item(i, 3);
        nameItem->setText(displayNameFull);
        nameItem->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter); // 完全居中
        nameItem->setFlags(nameItem->flags() & ~Qt::ItemIsEditable);
        nameItem->setData(Qt::UserRole, p.name);

        // --- 第 4 列: 下限 (调节框) ---
        SmartDoubleSpinBox* spinMin = new SmartDoubleSpinBox();
        spinMin->setStyleSheet(transparentSpinStyle); // 完全透明不遮挡
        spinMin->setRange(-9e9, 9e9);
        spinMin->setDecimals(10);
        spinMin->setValue(p.min);
        spinMin->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 4, spinMin);

        // --- 第 5 列: 上限 (调节框) ---
        SmartDoubleSpinBox* spinMax = new SmartDoubleSpinBox();
        spinMax->setStyleSheet(transparentSpinStyle); // 完全透明不遮挡
        spinMax->setRange(-9e9, 9e9);
        spinMax->setDecimals(10);
        spinMax->setValue(p.max);
        spinMax->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 5, spinMax);

        // --- 第 6 列: 步长 (调节框) ---
        SmartDoubleSpinBox* spinStep = new SmartDoubleSpinBox();
        spinStep->setStyleSheet(transparentSpinStyle); // 完全透明不遮挡
        spinStep->setRange(0.0, 10000.0);
        spinStep->setDecimals(10);
        spinStep->setValue(p.step);
        spinStep->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 6, spinStep);

        // --- 第 7 列: 显示 (复选框) ---
        QWidget* pWidgetVis = new QWidget();
        pWidgetVis->setStyleSheet(transparentWidgetStyle);
        QHBoxLayout* pLayoutVis = new QHBoxLayout(pWidgetVis);
        QCheckBox* chkVis = new QCheckBox();
        chkVis->setChecked(p.isVisible);
        chkVis->setStyleSheet(checkBoxStyle);
        pLayoutVis->addWidget(chkVis);
        pLayoutVis->setAlignment(Qt::AlignCenter);
        pLayoutVis->setContentsMargins(0,0,0,0);
        ui->tableWidget->setCellWidget(i, 7, pWidgetVis);

        // 联动逻辑：当勾选“拟合变量”(chkFit) 时，强制显示(chkVis)，并调用更新外观
        connect(chkFit, &QCheckBox::checkStateChanged, this, [this, chkVis, i](Qt::CheckState state){
            bool isFit = (state == Qt::Checked);
            if (isFit) {
                chkVis->setChecked(true);
                chkVis->setEnabled(false);
                chkVis->setStyleSheet("QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #ccc; border-radius: 3px; background-color: #e0e0e0; } "
                                      "QCheckBox::indicator:checked { background-color: #80bbeb; border-color: #80bbeb; }");
            } else {
                chkVis->setEnabled(true);
                chkVis->setStyleSheet("QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #cccccc; border-radius: 3px; background-color: white; }"
                                      "QCheckBox::indicator:checked { background-color: #0078d7; border-color: #0078d7; }"
                                      "QCheckBox::indicator:hover { border-color: #0078d7; }");
            }
            this->updateRowAppearance(i, isFit);
        });

        // 初始化关联状态
        if (p.isFit) {
            chkVis->setChecked(true);
            chkVis->setEnabled(false);
            chkVis->setStyleSheet("QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #ccc; border-radius: 3px; background-color: #e0e0e0; } "
                                  "QCheckBox::indicator:checked { background-color: #80bbeb; border-color: #80bbeb; }");
        }

        // 初始化加载时渲染该行的背景外观
        updateRowAppearance(i, p.isFit);
    }

    ui->tableWidget->resizeColumnsToContents();
    ui->tableWidget->horizontalHeader()->setSectionResizeMode(3, QHeaderView::Stretch); // 参数名称列自适应拉伸
}

// 内部函数：提取界面数据 (必须按新调整的列索引抓取数据)
void ParamSelectDialog::collectData()
{
    for(int i = 0; i < ui->tableWidget->rowCount(); ++i) {
        if(i >= m_params.size()) break;

        // 0列: 拟合变量
        QWidget* wFit = ui->tableWidget->cellWidget(i, 0);
        if (wFit) { QCheckBox* cb = wFit->findChild<QCheckBox*>(); if(cb) m_params[i].isFit = cb->isChecked(); }

        // 1列: 当前数值
        QDoubleSpinBox* spinVal = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 1));
        if(spinVal) m_params[i].value = spinVal->value();

        // 4列: 下限
        QDoubleSpinBox* spinMin = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 4));
        if(spinMin) m_params[i].min = spinMin->value();

        // 5列: 上限
        QDoubleSpinBox* spinMax = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 5));
        if(spinMax) m_params[i].max = spinMax->value();

        // 6列: 步长
        QDoubleSpinBox* spinStep = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 6));
        if(spinStep) m_params[i].step = spinStep->value();

        // 7列: 显示
        QWidget* wVis = ui->tableWidget->cellWidget(i, 7);
        if (wVis) { QCheckBox* cb = wVis->findChild<QCheckBox*>(); if(cb) m_params[i].isVisible = cb->isChecked(); }
    }
}

// 供外部调用
QList<FitParameter> ParamSelectDialog::getUpdatedParams() const
{
    return m_params;
}

double ParamSelectDialog::getFittingTime() const
{
    return ui->spinTimeMax->value();
}

void ParamSelectDialog::onConfirm()
{
    collectData();
    accept();
}

void ParamSelectDialog::onCancel()
{
    reject();
}

