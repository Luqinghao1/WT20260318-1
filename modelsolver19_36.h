/*
 * 文件名: modelsolver19_36.h
 * 文件作用与功能描述:
 * 1. 本代码文件为压裂水平井页岩型复合模型 Group 2 (Model 37-72) 的核心计算类头文件。
 * 2. 负责计算软件定义的另外 36 个以页岩型介质为主的试井模型理论图版数据。
 * 3. 涵盖三种储层组合的数学求解:
 * - 页岩型 + 页岩型   (Model 37-48 / 内部编号 1-12)
 * - 页岩型 + 均质     (Model 49-60 / 内部编号 13-24)
 * - 页岩型 + 双重孔隙 (Model 61-72 / 内部编号 25-36)
 * 4. 实现了用于描述页岩型介质的瞬态平板模型 (Transient Slab Matrix)。
 * 5. 基于与 Group 1 一致的高稳定性 BEM 算法，支持 Fair 和 Hegeman 模型的 Laplace 域变井储耦合。
 */

#ifndef MODELSOLVER19_36_H
#define MODELSOLVER19_36_H

#include <QMap>
#include <QVector>
#include <QString>
#include <tuple>
#include <functional>
#include <QtConcurrent>

// 类型定义: <时间序列(t), 压力序列(Dp), 导数序列(Dp')>
using ModelCurveData = std::tuple<QVector<double>, QVector<double>, QVector<double>>;

class ModelSolver19_36
{
public:
    // 模型类型枚举 (对应 modelsolver2.csv 中的 model2-1 到 model2-36)
    // 内部索引 0-35 (对应全局模型 ID 37-72)，保持数量和对应关系不变
    enum ModelType {
        // --- 第一组: 页岩型 + 页岩型 (37-48) ---
        // 无限大 (Const, Line, Fair, Hegeman)
        Model_1 = 0, Model_2, Model_3, Model_4,
        // 封闭
        Model_5, Model_6, Model_7, Model_8,
        // 定压
        Model_9, Model_10, Model_11, Model_12,

        // --- 第二组: 页岩型 + 均质 (49-60) ---
        Model_13, Model_14, Model_15, Model_16, // 无限大
        Model_17, Model_18, Model_19, Model_20, // 封闭
        Model_21, Model_22, Model_23, Model_24, // 定压

        // --- 第三组: 页岩型 + 双重孔隙 (61-72) ---
        Model_25, Model_26, Model_27, Model_28, // 无限大
        Model_29, Model_30, Model_31, Model_32, // 封闭
        Model_33, Model_34, Model_35, Model_36  // 定压
    };

    /**
     * @brief 构造函数，初始化特定类型的求解器
     * @param type 模型类型
     */
    explicit ModelSolver19_36(ModelType type);

    /**
     * @brief 析构函数，释放相关资源
     */
    virtual ~ModelSolver19_36();

    /**
     * @brief 设置高精度计算模式
     */
    void setHighPrecision(bool high);

    /**
     * @brief 计算理论曲线的核心接口
     * @param params 模型参数集合 (包含 k, S, cD, alpha, C_phi 等物理参数)
     * @param providedTime 指定的时间序列 (若为空则由内部自动生成对数时间分布)
     * @return 返回包含时间、压力、导数的元组
     */
    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime = QVector<double>());

    /**
     * @brief 获取模型的中文名称和描述，用于界面显示
     */
    static QString getModelName(ModelType type, bool verbose = true);

    /**
     * @brief 生成对数分布的时间步，满足早期至晚期的全阶段时间跨度
     */
    static QVector<double> generateLogTimeSteps(int count, double startExp, double endExp);

private:
    /**
     * @brief 驱动 Stehfest 反演算法，计算无因次压力和导数
     */
    void calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& params,
                             std::function<double(double, const QMap<QString, double>&)> laplaceFunc,
                             QVector<double>& outPD, QVector<double>& outDeriv);

    /**
     * @brief Laplace空间下的复合模型总控制与井储耦合计算函数
     */
    double flaplace_composite(double z, const QMap<QString, double>& p);

    /**
     * @brief 核心边界元方法(BEM)求解函数，构建并求解多段压裂水平井压力分布
     * 注意：去除了 spacingD，改为内部自动生成 [0.1, 0.9] 坐标系
     */
    double PWD_composite(double z, double fs1, double fs2, double M12, double LfD, double rmD, double reD,
                         int n_seg, int n_fracs, ModelType type);

    // --- 介质函数 ---
    /**
     * @brief 计算双重孔隙介质的拉普拉斯域函数 f(s)
     */
    double calc_fs_dual(double u, double omega, double lambda);

    /**
     * @brief 计算页岩(瞬态平板)介质的拉普拉斯域函数 f(s)
     */
    double calc_fs_shale(double u, double omega, double lambda);

    // --- 数学辅助 ---
    /**
     * @brief 标度修正的 Bessel I 函数 (消去大实参时的指数爆炸)
     */
    double scaled_besseli(int v, double x);

    /**
     * @brief 标度修正的 Bessel K 函数 (配合指数修正以杜绝溢出崩溃)
     */
    static double safe_bessel_k_scaled(int v, double x);

    double gauss15(std::function<double(double)> f, double a, double b);
    double adaptiveGauss(std::function<double(double)> f, double a, double b, double eps, int depth, int maxDepth);

    // --- Stehfest 反演算法辅助 ---
    double getStehfestCoeff(int i, int N);
    void precomputeStehfestCoeffs(int N);
    double factorial(int n);

private:
    ModelType m_type;
    bool m_highPrecision;
    QVector<double> m_stehfestCoeffs;
    int m_currentN;
};

#endif // MODELSOLVER19_36_H
