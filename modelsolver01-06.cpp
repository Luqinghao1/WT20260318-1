/*
 * 文件名: modelsolver01-06.cpp
 * 文件作用与功能描述:
 * 1. 压裂水平井复合模型 Group 1 (Model 1-36) 的核心计算算法实现。
 * 2. 依据最新的理论模型推导，将水平井坐标系原点设定为井筒起始点。
 * 3. 实现了高稳定性的边界元方法 (BEM)，采用标度修正的 Bessel 函数计算，彻底杜绝了大参数下的 NaN 或 Inf。
 * 4. [物理修正] 将变井储效应（Fair 模型、Hegeman 模型）从时间域移至 Laplace 域进行严谨耦合计算。
 * 5. [异常修复] 通过分离含有对数奇异性（log(0)）的对角积分部分解析处理，消除了积分杂音和微小波纹。
 * 6. 使用增大的 L_spacing (Bourdet 平滑跨度)，有效消除了压力导数的高频数值震荡。
 */

#include "modelsolver01-06.h"
#include "pressurederivativecalculator.h"

#include <Eigen/Dense>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <QDebug>
#include <QtConcurrent>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 二维点结构体，用于边界元离散化空间定位
struct Point2D { double x; double y; };

// 基础的 Bessel K 函数，增加下限托底防异常
static double safe_bessel_k(int v, double x) {
    if (x < 1e-15) x = 1e-15;
    try { return boost::math::cyl_bessel_k(v, x); } catch (...) { return 0.0; }
}

// 标度修正的 Bessel K 函数: 返回 K_v(x) * exp(x)，用于防止大实参下的数值溢出或直接返回0
double ModelSolver01_06::safe_bessel_k_scaled(int v, double x) {
    if (x < 1e-15) x = 1e-15;
    // 当 x 很大时 (> 600)，调用渐近展开以避免在未乘 exp(x) 之前产生双精度浮点数 underflow
    if (x > 600.0) return std::sqrt(M_PI / (2.0 * x));
    try { return boost::math::cyl_bessel_k(v, x) * std::exp(x); } catch (...) { return 0.0; }
}

// 标度修正的 Bessel I 函数: 返回 I_v(x) * exp(-x)，用于防范指数爆炸
static double safe_bessel_i_scaled(int v, double x) {
    if (x < 0) x = -x;
    if (x > 600.0) return 1.0 / std::sqrt(2.0 * M_PI * x); // 渐近展开近似
    try { return boost::math::cyl_bessel_i(v, x) * std::exp(-x); } catch (...) { return 0.0; }
}

// 构造函数，初始化 Stehfest 数组
ModelSolver01_06::ModelSolver01_06(ModelType type)
    : m_type(type), m_highPrecision(true), m_currentN(0) {
    // 提升 Stehfest 阶数至 8 以增强复相介质和变井储条件下的抗震荡性
    precomputeStehfestCoeffs(8);
}

ModelSolver01_06::~ModelSolver01_06() {}

// 启用高精度计算选项设定
void ModelSolver01_06::setHighPrecision(bool high) { m_highPrecision = high; }

// 根据模型编号解析对应的界面显示名称
QString ModelSolver01_06::getModelName(ModelType type, bool verbose)
{
    int id = (int)type + 1; // 转换为 1-36 的 ID
    QString baseName;
    QString subType;

    // 确定储层类型 (每12个一组)
    if (id <= 12) {
        baseName = QString("夹层型储层试井解释模型%1").arg(id);
        subType = "夹层型+夹层型";
    } else if (id <= 24) {
        baseName = QString("夹层型储层试井解释模型%1").arg(id);
        subType = "夹层型+均质";
    } else {
        baseName = QString("径向复合模型%1").arg(id - 24); // 径向复合从1开始编号
        subType = "均质+均质";
    }

    if (!verbose) return baseName;

    // 解析井筒储集模型分类 (0: 定井储, 1: 线源解, 2: Fair, 3: Hegeman)
    int rem4 = (id - 1) % 4;
    QString strStorage;
    if (rem4 == 0) strStorage = "定井储";
    else if (rem4 == 1) strStorage = "线源解";
    else if (rem4 == 2) strStorage = "Fair模型";
    else strStorage = "Hegeman模型";

    // 解析外边界条件类别
    int groupIdx = (id - 1) % 12;
    QString strBoundary;
    if (groupIdx < 4) strBoundary = "无限大外边界";
    else if (groupIdx < 8) strBoundary = "封闭边界";
    else strBoundary = "定压边界";

    return QString("%1\n(%2、%3、%4)").arg(baseName).arg(strStorage).arg(strBoundary).arg(subType);
}

// 生成对数时间节点，提供用于反演的数据坐标
QVector<double> ModelSolver01_06::generateLogTimeSteps(int count, double startExp, double endExp)
{
    QVector<double> t;
    if (count <= 0) return t;
    t.reserve(count);
    for (int i = 0; i < count; ++i) {
        double exponent = startExp + (endExp - startExp) * i / (count - 1);
        t.append(std::pow(10.0, exponent));
    }
    return t;
}

// 理论压力与压力导数计算的主控流程
ModelCurveData ModelSolver01_06::calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime)
{
    QVector<double> tPoints = providedTime;
    if (tPoints.isEmpty()) tPoints = generateLogTimeSteps(100, -3.0, 4.0);

    // --- 1. 读取基础物理参数 ---
    double phi = params.value("phi", 0.05);
    double mu = params.value("mu", 0.5);
    double B = params.value("B", 1.2);
    double Ct = params.value("Ct", 5e-4);
    double q = params.value("q", 50.0);
    double h = params.value("h", 20.0);
    double kf = params.value("kf", 1e-2);
    double L = params.value("L", 1000.0);

    if (L < 1e-9) L = 1000.0;
    // 参数有效性拦截，防止除零
    if (phi < 1e-12 || mu < 1e-12 || Ct < 1e-12 || kf < 1e-12) {
        return std::make_tuple(tPoints, QVector<double>(tPoints.size(), 0.0), QVector<double>(tPoints.size(), 0.0));
    }

    // --- 2. 转换无因次时间序列 ---
    double td_coeff = 14.4 * kf / (phi * mu * Ct * std::pow(L, 2.0));
    QVector<double> tD_vec;
    tD_vec.reserve(tPoints.size());
    for(double t : tPoints) tD_vec.append(td_coeff * t);

    // --- 3. Stehfest 参数设定 ---
    QMap<QString, double> calcParams = params;
    int N = (int)calcParams.value("N", 8); // 强制匹配优化的 8 阶数
    if (N < 4 || N > 18 || N % 2 != 0) N = 8;
    calcParams["N"] = N;
    precomputeStehfestCoeffs(N);

    if (!calcParams.contains("nf") || calcParams["nf"] < 1) calcParams["nf"] = 4;
    if (!calcParams.contains("n_seg")) calcParams["n_seg"] = 6;

    // --- 4. Stehfest 数值反演调用 (将变井储留到 Laplace 中解算) ---
    QVector<double> PD_vec, Deriv_vec;
    auto func = std::bind(&ModelSolver01_06::flaplace_composite, this, std::placeholders::_1, std::placeholders::_2);
    calculatePDandDeriv(tD_vec, calcParams, func, PD_vec, Deriv_vec);

    // --- 5. 获取有因次压力结果 ---
    double p_coeff = 1.842e-3 * q * mu * B / (kf * h);
    QVector<double> finalP(tPoints.size()), finalDP(tPoints.size());

    for(int i = 0; i < tPoints.size(); ++i) {
        finalP[i] = p_coeff * PD_vec[i];
    }

    // --- 6. 计算平滑压力导数 (使用扩大的跨度消除震荡) ---
    if (tPoints.size() > 2) {
        // [调整] L_spacing 修改为 0.4，以彻底过滤高频数值波纹
        finalDP = PressureDerivativeCalculator::calculateBourdetDerivative(tPoints, finalP, 0.4);
    } else {
        finalDP.fill(0.0);
    }

    return std::make_tuple(tPoints, finalP, finalDP);
}

// 基于 Stehfest 的时间域数值反演具体实现
void ModelSolver01_06::calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& params,
                                           std::function<double(double, const QMap<QString, double>&)> laplaceFunc,
                                           QVector<double>& outPD, QVector<double>& outDeriv)
{
    int numPoints = tD.size();
    outPD.resize(numPoints);
    outDeriv.resize(numPoints);

    int N = (int)params.value("N", 8);
    double ln2 = 0.6931471805599453;
    double gamaD = params.value("gamaD", 0.002); // 压敏系数提取

    QVector<int> indexes(numPoints);
    std::iota(indexes.begin(), indexes.end(), 0);

    // 单点无因次压力反演匿名函数
    auto calculateSinglePoint = [&](int k) {
        double t = tD[k];
        if (t <= 1e-10) { outPD[k] = 0.0; return; }

        double pd_val = 0.0;
        for (int m = 1; m <= N; ++m) {
            double z = m * ln2 / t;
            double pf = laplaceFunc(z, params); // 获取耦合井筒的 Laplace 压力
            if (std::isnan(pf) || std::isinf(pf)) pf = 0.0;
            pd_val += getStehfestCoeff(m, N) * pf;
        }

        double pd_real = pd_val * ln2 / t;

        // --- 摄动法压敏修正与防复数保护 ---
        if (gamaD > 1e-9) {
            double arg = 1.0 - gamaD * pd_real;
            // 严禁 arg 小于等于 0 导致对数抛出复数错误引发的图版倒V尖点
            if (arg <= 1e-6) arg = 1e-6;
            pd_real = -1.0 / gamaD * std::log(arg);
        }

        // 托底保护，确保无负数值导致双对数图不可视
        if (pd_real <= 1e-15) pd_real = 1e-15;
        outPD[k] = pd_real;
    };

    // 启用多线程并发加速计算
    QtConcurrent::blockingMap(indexes, calculateSinglePoint);
}

// Laplace 空间下耦合了边界元和井储效应的总函数
double ModelSolver01_06::flaplace_composite(double z, const QMap<QString, double>& p) {
    // 基础尺寸提取
    double M12 = p.value("M12", 10.0);
    double L = p.value("L", 800.0);
    double Lf = p.value("Lf", 40.0);
    double rm = p.value("rm", 1000.0);
    double re = p.value("re", 20000.0);
    double eta12 = p.value("eta12", 0.2);

    // 无因次化距离
    double LfD = (L > 1e-9) ? Lf / L : 0.05;
    double rmD = (L > 1e-9) ? rm / L : 1.25;
    double reD = (L > 1e-9) ? re / L : 25.0;

    int n_fracs = (int)p.value("nf", 4);
    int n_seg = (int)p.value("n_seg", 6);

    double fs1 = 1.0, fs2 = 1.0;
    int id = (int)m_type + 1;

    // --- 定义复相介质函数 f(s) ---
    // 组别 1 & 2: 内部双重介质
    if (id <= 24) {
        double omga1 = p.value("omega1", 0.4);
        double remda1 = p.value("lambda1", 1e-3);
        double one_minus = 1.0 - omga1;
        fs1 = (omga1 * one_minus * z + remda1) / (one_minus * z + remda1);
    }
    // 组别 1: 外部双重介质
    if (id <= 12) {
        double omga2 = p.value("omega2", 0.08);
        double remda2 = p.value("lambda2", 1e-4);
        double one_minus = 1.0 - omga2;
        fs2 = eta12 * (omga2 * one_minus * eta12 * z + remda2) / (one_minus * eta12 * z + remda2);
    } else {
        fs2 = eta12; // 均质外区
    }

    // --- 调用 BEM 核心获取纯地层解 ---
    double pf_base = PWD_composite(z, fs1, fs2, M12, LfD, rmD, reD, n_seg, n_fracs, m_type);

    // --- Laplace 域变井筒储集耦合 (杜哈美原理物理修正式) ---
    int storageType = (int)m_type % 4; // 0:定井储, 1:线源解, 2:Fair模型, 3:Hegeman模型

    // 若为线源解，直接返回且不加任何表皮与井储
    if (storageType == 1) {
        return pf_base;
    }

    // 计算包含 Fair / Hegeman 的等效动态井筒储集系数 (sCa)
    double CD = p.value("cD", 0.15); // 无因次井筒储存系数
    double S = p.value("S", 10.0);
    if (S < 0.0) S = 0.0;

    double alpha = p.value("alpha", 1e-1);
    double C_phi = p.value("C_phi", 1e-4);
    double sCa = CD;

    if (storageType == 2) {
        // Fair 变井储解析公式
        if (std::abs(alpha) > 1e-12) {
            sCa = CD * (1.0 + C_phi / (1.0 + alpha * z));
        }
    } else if (storageType == 3) {
        // Hegeman 变井储误差解析公式
        if (std::abs(alpha) > 1e-12) {
            double xx = z * alpha / 2.0;
            double exp_erfc = 0.0;
            // 防溢出保护：利用渐近线处理极大自变量下的 erfc 计算
            if (xx < 10.0) {
                exp_erfc = std::exp(xx * xx) * std::erfc(xx);
            } else {
                exp_erfc = 1.0 / (std::sqrt(M_PI) * xx);
            }
            sCa = CD * (1.0 + C_phi * exp_erfc);
        }
    }

    // 耦合计算公式合并
    double pf = 0.0;
    double num = z * pf_base + S;
    double den = z + sCa * z * z * num;
    if (std::abs(den) > 1e-100) {
        pf = num / den;
    } else {
        pf = pf_base; // 极端防御性降级
    }

    return pf;
}

// 基于 BEM 的水平井面源积分方程计算
double ModelSolver01_06::PWD_composite(double z, double fs1, double fs2, double M12, double LfD, double rmD, double reD,
                                       int n_seg, int n_fracs, ModelType type) {
    // 判定外边界条件类型
    int id = (int)type + 1;
    int groupIdx = (id - 1) % 12;
    bool isInfinite = (groupIdx < 4);
    bool isClosed = (groupIdx >= 4 && groupIdx < 8);
    bool isConstP = (groupIdx >= 8);

    double gama1 = std::sqrt(z * fs1);
    double gama2 = std::sqrt(z * fs2);
    double arg_g1_rm = gama1 * rmD;
    double arg_g2_rm = gama2 * rmD;

    // --- 采用缩放修正公式获取界面耦合贝塞尔系数，绝对阻断大自变量极化 ---
    double k0_g1 = safe_bessel_k_scaled(0, arg_g1_rm);
    double k1_g1 = safe_bessel_k_scaled(1, arg_g1_rm);
    double i0_g1 = safe_bessel_i_scaled(0, arg_g1_rm);
    double i1_g1 = safe_bessel_i_scaled(1, arg_g1_rm);

    double k0_g2 = safe_bessel_k_scaled(0, arg_g2_rm);
    double k1_g2 = safe_bessel_k_scaled(1, arg_g2_rm);
    double i0_g2 = safe_bessel_i_scaled(0, arg_g2_rm);
    double i1_g2 = safe_bessel_i_scaled(1, arg_g2_rm);

    double T1_prime = k0_g2;
    double T2_prime = -k1_g2;

    // 非无限大外边界反射干涉耦合
    if (!isInfinite && reD > 1e-5) {
        double arg_re = gama2 * reD;
        double k0_re = safe_bessel_k_scaled(0, arg_re);
        double k1_re = safe_bessel_k_scaled(1, arg_re);
        double i0_re = safe_bessel_i_scaled(0, arg_re);
        double i1_re = safe_bessel_i_scaled(1, arg_re);

        // 指数平移项：因为 rmD <= reD，其指数为负值可安全衰减
        double exp_factor = std::exp(2.0 * gama2 * (rmD - reD));

        if (isClosed) {
            double ratio = k1_re / std::max(i1_re, 1e-100);
            T1_prime = ratio * i0_g2 * exp_factor + k0_g2;
            T2_prime = ratio * i1_g2 * exp_factor - k1_g2;
        } else if (isConstP) {
            double ratio = -k0_re / std::max(i0_re, 1e-100);
            T1_prime = ratio * i0_g2 * exp_factor + k0_g2;
            T2_prime = ratio * i1_g2 * exp_factor - k1_g2;
        }
    }

    double Acup_prime   = M12 * gama1 * k1_g1 * T1_prime + gama2 * k0_g1 * T2_prime;
    double Acdown_prime = M12 * gama1 * i1_g1 * T1_prime - gama2 * i0_g1 * T2_prime;
    if (std::abs(Acdown_prime) < 1e-100) {
        Acdown_prime = (Acdown_prime >= 0) ? 1e-100 : -1e-100;
    }
    double Ac_core = Acup_prime / Acdown_prime;

    // --- 横向多段压裂坐标系重塑 (以起始点为原点) ---
    int total_segments = n_fracs * n_seg;
    double segLen = 2.0 * LfD / n_seg;
    double halfLen = segLen / 2.0;

    QVector<Point2D> segmentCenters;
    segmentCenters.reserve(total_segments);

    // 生成横贯整个水平井长度 L 区间的分布位置 (0.1, ..., 0.9)
    QVector<double> xwD(n_fracs);
    if (n_fracs == 1) {
        xwD[0] = 0.5;
    } else {
        for (int k = 0; k < n_fracs; ++k) {
            xwD[k] = 0.1 + 0.8 * (double)k / (double)(n_fracs - 1);
        }
    }

    // 挂载所有裂缝子微元的二维平面坐标
    for (int i = 0; i < n_fracs; ++i) {
        for (int j = 0; j < n_seg; ++j) {
            double currentX = xwD[i];
            double currentY = -LfD + (j + 0.5) * segLen;
            segmentCenters.append({currentX, currentY});
        }
    }

    // --- BEM 矩阵注入 ---
    int size = total_segments + 1;
    Eigen::MatrixXd A_mat(size, size);
    Eigen::VectorXd b_vec(size);
    b_vec.setZero();

    for (int i = 0; i < total_segments; ++i) {
        for (int j = i; j < total_segments; ++j) {
            Point2D pi = segmentCenters[i];
            Point2D pj = segmentCenters[j];
            double dx_sq = (pi.x - pj.x) * (pi.x - pj.x);

            double val = 0.0;
            if (i == j) {
                // [优化处理] 对角积分：提取出对数奇点以做精确的解析累加，避开数值波纹陷阱
                auto diag_integrand = [&](double a) -> double {
                    if (a < 1e-14) {
                        return -std::log(gama1 / 2.0) - 0.57721566490153286 + Ac_core * safe_bessel_i_scaled(0, 0) * std::exp(-2.0 * arg_g1_rm);
                    }
                    double arg_dist = gama1 * a;
                    // 原式：K_0 + 解析提走的 log(a) + I_0 复合反射
                    return safe_bessel_k(0, arg_dist) + std::log(a) + Ac_core * safe_bessel_i_scaled(0, arg_dist) * std::exp(arg_dist - 2.0 * arg_g1_rm);
                };

                double val_num = adaptiveGauss(diag_integrand, 0.0, halfLen, 1e-10, 0, 15);
                double val_ana = halfLen * (1.0 - std::log(halfLen));
                val = 2.0 * (val_num + val_ana);
            } else {
                // 异位置的非奇异积分
                auto offdiag_integrand = [&](double a) -> double {
                    double dy = pi.y - (pj.y + a);
                    double dist_val = std::sqrt(dx_sq + dy * dy);
                    double arg_dist = gama1 * dist_val;

                    double term1 = safe_bessel_k(0, arg_dist);
                    double exponent = arg_dist - 2.0 * arg_g1_rm;
                    double term2 = Ac_core * safe_bessel_i_scaled(0, arg_dist) * std::exp(exponent);
                    return term1 + term2;
                };
                val = adaptiveGauss(offdiag_integrand, -halfLen, halfLen, 1e-10, 0, 15);
            }

            double element = val / (M12 * segLen);
            A_mat(i, j) = element;
            if (i != j) A_mat(j, i) = element;
        }
    }

    // 矩阵末端强制条件引入与归一化缩放 (阻挡病态方程条件数爆炸)
    for (int i = 0; i < total_segments; ++i) {
        A_mat(i, total_segments) = -1.0;
        A_mat(total_segments, i) = 1.0;
    }
    A_mat(total_segments, total_segments) = 0.0;
    b_vec(total_segments) = 1.0 / z; // 保持等号两侧同步除以 z 进行数量级抹平

    // 使用全主元 LU 矩阵求解以防止矩阵结构性失真
    Eigen::VectorXd x_sol = A_mat.fullPivLu().solve(b_vec);

    // 提取计算得到的复合系统主控压力节点
    return x_sol(total_segments);
}

// --- 基础数学子程序实现 ---
double ModelSolver01_06::scaled_besseli(int v, double x) { return safe_bessel_i_scaled(v, x); }

// Gauss-Kronrod 15 点定积分基础单元
double ModelSolver01_06::gauss15(std::function<double(double)> f, double a, double b) {
    static const double X[] = { 0.0, 0.20119409, 0.39415135, 0.57097217, 0.72441773, 0.84820658, 0.93729853, 0.98799252 };
    static const double W[] = { 0.20257824, 0.19843149, 0.18616100, 0.16626921, 0.13957068, 0.10715922, 0.07036605, 0.03075324 };
    double h = 0.5 * (b - a);
    double c = 0.5 * (a + b);
    double s = W[0] * f(c);
    for (int i = 1; i < 8; ++i) {
        double dx = h * X[i];
        s += W[i] * (f(c - dx) + f(c + dx));
    }
    return s * h;
}

// 递进式自适应递归积分包装器
double ModelSolver01_06::adaptiveGauss(std::function<double(double)> f, double a, double b, double eps, int depth, int maxDepth) {
    double c = (a + b) / 2.0;
    double v1 = gauss15(f, a, b);
    double v2 = gauss15(f, a, c) + gauss15(f, c, b);
    if (depth >= maxDepth || std::abs(v1 - v2) < eps * (std::abs(v2) + 1.0)) return v2;
    return adaptiveGauss(f, a, c, eps/2, depth+1, maxDepth) + adaptiveGauss(f, c, b, eps/2, depth+1, maxDepth);
}

// Stehfest 累加系数池预运算加载
void ModelSolver01_06::precomputeStehfestCoeffs(int N) {
    if (m_currentN == N && !m_stehfestCoeffs.isEmpty()) return;
    m_currentN = N; m_stehfestCoeffs.resize(N + 1);
    for (int i = 1; i <= N; ++i) {
        double s = 0.0;
        int k1 = (i + 1) / 2;
        int k2 = std::min(i, N / 2);
        for (int k = k1; k <= k2; ++k) {
            double num = std::pow((double)k, N / 2.0) * factorial(2 * k);
            double den = factorial(N / 2 - k) * factorial(k) * factorial(k - 1) * factorial(i - k) * factorial(2 * k - i);
            if (den != 0) s += num / den;
        }
        double sign = ((i + N / 2) % 2 == 0) ? 1.0 : -1.0;
        m_stehfestCoeffs[i] = sign * s;
    }
}

// 安全提取反演常数
double ModelSolver01_06::getStehfestCoeff(int i, int N) {
    if (m_currentN != N) return 0.0;
    if (i < 1 || i > N) return 0.0;
    return m_stehfestCoeffs[i];
}

// 高速阶乘运算
double ModelSolver01_06::factorial(int n) {
    if(n <= 1) return 1.0;
    double r = 1.0;
    for(int i = 2; i <= n; ++i) r *= i;
    return r;
}
