/*
 * 文件名: fittingcore.cpp
 * 文件作用: 试井拟合核心算法实现
 */

#include "fittingcore.h"
#include "modelparameter.h" // 引入模型参数单例
#include <QtConcurrent>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <Eigen/Dense>

FittingCore::FittingCore(QObject *parent)
    : QObject(parent), m_modelManager(nullptr), m_isCustomSamplingEnabled(false), m_stopRequested(false)
{
    // 监听异步任务完成
    connect(&m_watcher, &QFutureWatcher<void>::finished, this, &FittingCore::sigFitFinished);
}

void FittingCore::setModelManager(ModelManager *m) {
    m_modelManager = m;
}

void FittingCore::setObservedData(const QVector<double> &t, const QVector<double> &p, const QVector<double> &d) {
    m_obsTime = t;
    m_obsDeltaP = p;
    m_obsDerivative = d;
}

void FittingCore::setSamplingSettings(const QList<SamplingInterval> &intervals, bool enabled) {
    m_customIntervals = intervals;
    m_isCustomSamplingEnabled = enabled;
}

void FittingCore::startFit(ModelManager::ModelType modelType, const QList<FitParameter> &params, double weight) {
    if (m_watcher.isRunning()) return;

    m_stopRequested = false;
    // 启动异步线程执行拟合
    m_watcher.setFuture(QtConcurrent::run([this, modelType, params, weight](){
        runOptimizationTask(modelType, params, weight);
    }));
}

void FittingCore::stopFit() {
    m_stopRequested = true;
}

void FittingCore::getLogSampledData(const QVector<double>& srcT, const QVector<double>& srcP, const QVector<double>& srcD,
                                    QVector<double>& outT, QVector<double>& outP, QVector<double>& outD)
{
    outT.clear(); outP.clear(); outD.clear();
    if (srcT.isEmpty()) return;

    // 辅助结构体用于排序去重
    struct DataPoint {
        double t, p, d;
        bool operator<(const DataPoint& other) const { return t < other.t; }
        bool operator==(const DataPoint& other) const { return std::abs(t - other.t) < 1e-9; }
    };
    QVector<DataPoint> points;

    // 模式1：默认策略
    if (!m_isCustomSamplingEnabled) {
        int targetCount = 200;
        if (srcT.size() <= targetCount) {
            outT = srcT; outP = srcP; outD = srcD;
            return;
        }

        double tMin = srcT.first() <= 1e-10 ? 1e-4 : srcT.first();
        double tMax = srcT.last();
        double logMin = log10(tMin);
        double logMax = log10(tMax);
        double step = (logMax - logMin) / (targetCount - 1);

        int currentIndex = 0;
        for (int i = 0; i < targetCount; ++i) {
            double targetT = pow(10, logMin + i * step);
            double minDiff = 1e30;
            int bestIdx = currentIndex;
            while (currentIndex < srcT.size()) {
                double diff = std::abs(srcT[currentIndex] - targetT);
                if (diff < minDiff) { minDiff = diff; bestIdx = currentIndex; }
                else break;
                currentIndex++;
            }
            currentIndex = bestIdx;
            points.append({srcT[bestIdx],
                           (bestIdx<srcP.size()?srcP[bestIdx]:0.0),
                           (bestIdx<srcD.size()?srcD[bestIdx]:0.0)});
        }
    }
    // 模式2：自定义区间策略
    else {
        if (m_customIntervals.isEmpty()) {
            outT = srcT; outP = srcP; outD = srcD;
            return;
        }
        for (const auto& interval : m_customIntervals) {
            double tStart = interval.tStart;
            double tEnd = interval.tEnd;
            int count = interval.count;
            if (count <= 0) continue;

            auto itStart = std::lower_bound(srcT.begin(), srcT.end(), tStart);
            auto itEnd = std::upper_bound(srcT.begin(), srcT.end(), tEnd);
            int idxStart = std::distance(srcT.begin(), itStart);
            int idxEnd = std::distance(srcT.begin(), itEnd);

            if (idxStart >= srcT.size() || idxStart >= idxEnd) continue;

            double subMin = srcT[idxStart];
            double subMax = srcT[idxEnd - 1];
            if (subMin <= 1e-10) subMin = 1e-4;

            double logMin = log10(subMin);
            double logMax = log10(subMax);
            double step = (count > 1) ? (logMax - logMin) / (count - 1) : 0;

            int subCurrentIdx = idxStart;
            for (int i = 0; i < count; ++i) {
                double targetT = (count == 1) ? subMin : pow(10, logMin + i * step);
                double minDiff = 1e30;
                int bestIdx = subCurrentIdx;
                while (subCurrentIdx < idxEnd) {
                    double diff = std::abs(srcT[subCurrentIdx] - targetT);
                    if (diff < minDiff) { minDiff = diff; bestIdx = subCurrentIdx; }
                    else break;
                    subCurrentIdx++;
                }
                subCurrentIdx = bestIdx;
                if (bestIdx < srcT.size()) {
                    points.append({srcT[bestIdx],
                                   (bestIdx<srcP.size()?srcP[bestIdx]:0.0),
                                   (bestIdx<srcD.size()?srcD[bestIdx]:0.0)});
                }
            }
        }
    }

    std::sort(points.begin(), points.end());
    auto last = std::unique(points.begin(), points.end());
    points.erase(last, points.end());

    for (const auto& p : points) {
        outT.append(p.t);
        outP.append(p.p);
        outD.append(p.d);
    }
}

// [核心修正] 参数预处理函数
QMap<QString, double> FittingCore::preprocessParams(const QMap<QString, double>& rawParams, ModelManager::ModelType type)
{
    QMap<QString, double> processed = rawParams;
    ModelParameter* mp = ModelParameter::instance();

    // 辅助 lambda：如果值不存在或接近0，则使用安全默认值
    auto getSafeParam = [&](const QString& key, double mpVal, double defaultVal) {
        if (rawParams.contains(key)) return rawParams[key];
        if (std::abs(mpVal) > 1e-15) return mpVal;
        return defaultVal;
    };

    // 1. 基础参数补全与校验
    double phi = getSafeParam("phi", mp->getPhi(), 0.05);
    double h   = getSafeParam("h",   mp->getH(),   20.0);
    double Ct  = getSafeParam("Ct",  mp->getCt(),  5e-4);
    double mu  = getSafeParam("mu",  mp->getMu(),  0.5);
    double B   = getSafeParam("B",   mp->getB(),   1.05);
    double q   = getSafeParam("q",   mp->getQ(),   5.0);
    double rw  = getSafeParam("rw",  mp->getRw(),  0.1);

    // 回写安全值
    processed["phi"] = phi;
    processed["h"] = h;
    processed["Ct"] = Ct;
    processed["mu"] = mu;
    processed["B"] = B;
    processed["q"] = q;
    processed["rw"] = rw;

    // 2. 处理 L 和 LfD
    double L = processed.value("L", 0.0);
    if (L < 1e-9) {
        L = 1000.0; // 默认 L
        processed["L"] = L;
    }
    // 计算 LfD
    if (processed.contains("Lf")) {
        processed["LfD"] = processed["Lf"] / L;
    } else {
        processed["LfD"] = 0.0;
    }

    // 3. 处理流度比 km -> M12
    if (!processed.contains("M12") && processed.contains("km")) {
        processed["M12"] = processed["km"];
    }

    // 4. 处理井筒储集 C -> cD
    bool hasStorage = (type == ModelManager::Model_1 || type == ModelManager::Model_3 || type == ModelManager::Model_5);
    if (hasStorage) {
        // 优先检查有因次 C，转换为 cD
        if (processed.contains("C")) {
            double valC = processed["C"];
            // [公式修正] CD = 0.159 * C / (phi * h * Ct * L^2)
            // 原代码使用 rw^2，现改为 L^2 (水平井长度)
            double denom = phi * h * Ct * L * L;
            double cD = 0.0;
            // 0.159 是 1/(2*pi) 的近似值
            if (denom > 1e-20) cD = 0.159 * valC / denom;
            processed["cD"] = cD;
        }
    } else {
        processed["cD"] = 0.0;
        processed["S"] = 0.0;
    }

    // 5. 补充边界参数 re
    bool isInfinite = (type == ModelManager::Model_1 || type == ModelManager::Model_2);
    if (isInfinite) {
        if (!processed.contains("re")) processed["re"] = 20000.0;
    }

    return processed;
}

void FittingCore::runOptimizationTask(ModelManager::ModelType modelType, QList<FitParameter> fitParams, double weight) {
    runLevenbergMarquardtOptimization(modelType, fitParams, weight);
}

void FittingCore::runLevenbergMarquardtOptimization(ModelManager::ModelType modelType, QList<FitParameter> params, double weight) {
    if(m_modelManager) m_modelManager->setHighPrecision(false);

    QVector<int> fitIndices;
    for(int i=0; i<params.size(); ++i) {
        if(params[i].isFit && params[i].name != "LfD") fitIndices.append(i);
    }
    int nParams = fitIndices.size();

    QMap<QString, double> currentParamMap;
    for(const auto& p : params) currentParamMap.insert(p.name, p.value);

    // [关键] 使用预处理后的参数进行初始计算
    QMap<QString, double> solverParams = preprocessParams(currentParamMap, modelType);

    QVector<double> fitT, fitP, fitD;
    getLogSampledData(m_obsTime, m_obsDeltaP, m_obsDerivative, fitT, fitP, fitD);

    QVector<double> residuals = calculateResiduals(currentParamMap, modelType, weight, fitT, fitP, fitD);
    double currentSSE = calculateSumSquaredError(residuals);

    ModelCurveData curve = m_modelManager->calculateTheoreticalCurve(modelType, solverParams);
    emit sigIterationUpdated(currentSSE/residuals.size(), currentParamMap, std::get<0>(curve), std::get<1>(curve), std::get<2>(curve));

    if(nParams == 0) {
        emit sigFitFinished();
        return;
    }

    double lambda = 0.01;
    int maxIter = 50;

    for(int iter = 0; iter < maxIter; ++iter) {
        if(m_stopRequested) break;
        if (!residuals.isEmpty() && (currentSSE / residuals.size()) < 3e-3) break;

        emit sigProgress(iter * 100 / maxIter);

        QVector<QVector<double>> J = computeJacobian(currentParamMap, residuals, fitIndices, modelType, params, weight, fitT, fitP, fitD);
        int nRes = residuals.size();

        QVector<QVector<double>> H(nParams, QVector<double>(nParams, 0.0));
        QVector<double> g(nParams, 0.0);

        for(int k=0; k<nRes; ++k) {
            for(int i=0; i<nParams; ++i) {
                g[i] += J[k][i] * residuals[k];
                for(int j=0; j<=i; ++j) {
                    H[i][j] += J[k][i] * J[k][j];
                }
            }
        }
        for(int i=0; i<nParams; ++i) {
            for(int j=i+1; j<nParams; ++j) H[i][j] = H[j][i];
        }

        bool stepAccepted = false;
        for(int tryIter=0; tryIter<5; ++tryIter) {
            QVector<QVector<double>> H_lm = H;
            for(int i=0; i<nParams; ++i) H_lm[i][i] += lambda * (1.0 + std::abs(H[i][i]));

            QVector<double> negG(nParams);
            for(int i=0;i<nParams;++i) negG[i] = -g[i];

            QVector<double> delta = solveLinearSystem(H_lm, negG);
            QMap<QString, double> trialMap = currentParamMap;

            for(int i=0; i<nParams; ++i) {
                int pIdx = fitIndices[i];
                QString pName = params[pIdx].name;
                double oldVal = currentParamMap[pName];
                bool isLog = (oldVal > 1e-12 && pName != "S" && pName != "nf");
                double newVal;
                if(isLog) newVal = pow(10.0, log10(oldVal) + delta[i]);
                else newVal = oldVal + delta[i];
                newVal = qMax(params[pIdx].min, qMin(newVal, params[pIdx].max));
                trialMap[pName] = newVal;
            }

            // [关键] 内部会自动调用 preprocessParams
            QVector<double> newRes = calculateResiduals(trialMap, modelType, weight, fitT, fitP, fitD);
            double newSSE = calculateSumSquaredError(newRes);

            if(newSSE < currentSSE) {
                currentSSE = newSSE;
                currentParamMap = trialMap;
                residuals = newRes;
                lambda /= 10.0;
                stepAccepted = true;

                // 更新曲线
                QMap<QString, double> trialSolverParams = preprocessParams(trialMap, modelType);
                ModelCurveData iterCurve = m_modelManager->calculateTheoreticalCurve(modelType, trialSolverParams);
                emit sigIterationUpdated(currentSSE/nRes, currentParamMap, std::get<0>(iterCurve), std::get<1>(iterCurve), std::get<2>(iterCurve));
                break;
            } else {
                lambda *= 10.0;
            }
        }
        if(!stepAccepted && lambda > 1e10) break;
    }

    if(m_modelManager) m_modelManager->setHighPrecision(true);

    QMap<QString, double> finalSolverParams = preprocessParams(currentParamMap, modelType);
    ModelCurveData finalCurve = m_modelManager->calculateTheoreticalCurve(modelType, finalSolverParams);
    emit sigIterationUpdated(currentSSE/residuals.size(), currentParamMap, std::get<0>(finalCurve), std::get<1>(finalCurve), std::get<2>(finalCurve));
}

QVector<double> FittingCore::calculateResiduals(const QMap<QString, double>& params, ModelManager::ModelType modelType, double weight,
                                                const QVector<double>& t, const QVector<double>& obsP, const QVector<double>& obsD) {
    if(!m_modelManager || t.isEmpty()) return QVector<double>();

    // [关键] 参数预处理
    QMap<QString, double> solverParams = preprocessParams(params, modelType);

    ModelCurveData res = m_modelManager->calculateTheoreticalCurve(modelType, solverParams, t);
    const QVector<double>& pCal = std::get<1>(res);
    const QVector<double>& dpCal = std::get<2>(res);

    QVector<double> r;
    double wp = weight;
    double wd = 1.0 - weight;

    int count = qMin((int)obsP.size(), (int)pCal.size());
    for(int i=0; i<count; ++i) {
        if(obsP[i] > 1e-10 && pCal[i] > 1e-10)
            r.append( (log(obsP[i]) - log(pCal[i])) * wp );
        else
            r.append(0.0);
    }
    int dCount = qMin((int)obsD.size(), (int)dpCal.size());
    dCount = qMin(dCount, count);
    for(int i=0; i<dCount; ++i) {
        if(obsD[i] > 1e-10 && dpCal[i] > 1e-10)
            r.append( (log(obsD[i]) - log(dpCal[i])) * wd );
        else
            r.append(0.0);
    }
    return r;
}

QVector<QVector<double>> FittingCore::computeJacobian(const QMap<QString, double>& params, const QVector<double>& baseResiduals,
                                                      const QVector<int>& fitIndices, ModelManager::ModelType modelType,
                                                      const QList<FitParameter>& currentFitParams, double weight,
                                                      const QVector<double>& t, const QVector<double>& obsP, const QVector<double>& obsD) {
    int nRes = baseResiduals.size();
    int nParams = fitIndices.size();
    QVector<QVector<double>> J(nRes, QVector<double>(nParams));

    QVector<int> indices(nParams);
    std::iota(indices.begin(), indices.end(), 0);

    auto computeColumn = [&](int j) -> QVector<double> {
        int idx = fitIndices[j];
        QString pName = currentFitParams[idx].name;
        double val = params.value(pName);
        bool isLog = (val > 1e-12 && pName != "S" && pName != "nf");

        double h;
        QMap<QString, double> pPlus = params;
        QMap<QString, double> pMinus = params;

        if(isLog) {
            h = 0.01;
            double valLog = log10(val);
            pPlus[pName] = pow(10.0, valLog + h);
            pMinus[pName] = pow(10.0, valLog - h);
        } else {
            h = 1e-4;
            pPlus[pName] = val + h;
            pMinus[pName] = val - h;
        }

        // calculateResiduals 内部会自动调用 preprocessParams，直接传递扰动参数即可
        QVector<double> rPlus = this->calculateResiduals(pPlus, modelType, weight, t, obsP, obsD);
        QVector<double> rMinus = this->calculateResiduals(pMinus, modelType, weight, t, obsP, obsD);

        QVector<double> col(nRes, 0.0);
        if(rPlus.size() == nRes && rMinus.size() == nRes) {
            for(int i=0; i<nRes; ++i) col[i] = (rPlus[i] - rMinus[i]) / (2.0 * h);
        }
        return col;
    };

    QList<QVector<double>> results = QtConcurrent::blockingMapped(indices, computeColumn);
    for(int j=0; j<nParams; ++j) {
        const QVector<double>& col = results[j];
        for(int i=0; i<nRes; ++i) {
            if (i < col.size()) J[i][j] = col[i];
        }
    }
    return J;
}

QVector<double> FittingCore::solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b) {
    int n = b.size();
    if (n == 0) return QVector<double>();
    Eigen::MatrixXd matA(n, n);
    Eigen::VectorXd vecB(n);
    for (int i = 0; i < n; ++i) {
        vecB(i) = b[i];
        for (int j = 0; j < n; ++j) matA(i, j) = A[i][j];
    }
    Eigen::VectorXd x = matA.ldlt().solve(vecB);
    QVector<double> res(n);
    for (int i = 0; i < n; ++i) res[i] = x(i);
    return res;
}

double FittingCore::calculateSumSquaredError(const QVector<double>& residuals) {
    double sse = 0.0;
    for(double v : residuals) sse += v*v;
    return sse;
}
