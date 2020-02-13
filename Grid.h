#ifndef GRID_H
#define GRID_H

#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <iostream>

#include "Axis.h"
#include "Tools.h"

using std::vector;
using std::string;
using std::ofstream;
using std::ifstream;
using std::reference_wrapper;

class HistogramBase
{
public:
    HistogramBase();
    HistogramBase(const vector<Axis>& ax);
    bool isInGrid(const vector<double>& pos) const;
    bool index(const vector<double>& pos, vector<size_t>& idx) const;
    bool address(const vector<double>& pos, size_t& addr) const;
    vector<Axis> getAxes() const;
    vector<vector<double>> getTable() const;
    double getGridSize() const;
    size_t getDimension() const;
protected:
    size_t                  mNDim;
    vector<Axis>            mAxes;
    vector<size_t>          mAccu;
    size_t                  mGridSize;
    vector<vector<double>>  mPointTable;

    void fillTable();
};

HistogramBase::HistogramBase():
mNDim(0), mAxes(0), mAccu(0), mGridSize(0), mPointTable(0) {}

HistogramBase::HistogramBase(const vector<Axis>& ax):
mNDim(ax.size()), mAxes(ax), mAccu(mNDim) {
    mAccu[0] = 1;
    mGridSize = 1;
    for (size_t i = 0; i < mNDim; ++i) {
        mAccu[i] = (i == 0) ? 1 : (mAccu[i - 1] * mAxes[i - 1].bin());
        mGridSize *= mAxes[i].bin();
    }
    fillTable();
}

void HistogramBase::fillTable() {
    vector<vector<double>> middlePoint(mNDim);
    for (size_t i = 0; i < mNDim; ++i) {
        middlePoint[i] = mAxes[i].middlePoint();
    }
    mPointTable.assign(mNDim, vector<double>(mGridSize, 0.0));
    for (size_t i = 0; i < mNDim; ++i) {
        size_t repeatAll = 1, repeatOne = 1;
        for (size_t j = i + 1; j < mNDim; ++j) {
            repeatOne *= middlePoint[j].size();
        }
        for (size_t j = 0; j < i; ++j) {
            repeatAll *= middlePoint[j].size();
        }
        const size_t in_i_sz = middlePoint[i].size();
        for (size_t l = 0; l < in_i_sz; ++l) {
            std::fill_n(begin(mPointTable[i]) + l * repeatOne, repeatOne, middlePoint[i][l]);
        }
        for (size_t k = 0; k < repeatAll - 1; ++k) {
        std::copy_n(begin(mPointTable[i]), repeatOne * in_i_sz,
                    begin(mPointTable[i]) + repeatOne * in_i_sz * (k + 1));
        }
    }
}

bool HistogramBase::isInGrid(const vector<double>& pos) const {
    auto it_val = pos.cbegin();
    auto it_ax = mAxes.cbegin();
    while (it_ax != mAxes.cend()) {
        if (!(it_ax->isInBoundary((*it_val)))) return false;
    }
    return true;
}

bool HistogramBase::index(const vector<double>& pos, vector<size_t>& idx) const {
    idx.assign(mNDim, 0);
    auto it_val = pos.cbegin();
    auto it_ax = mAxes.cbegin();
    auto it_idx = idx.begin();
    while (it_idx != idx.end()) {
        if (!(it_ax->index((*it_val), (*it_idx)))) return false;
    }
    return true;
}

bool HistogramBase::address(const vector<double>& pos, size_t& addr) const {
    addr = 0;
    for (size_t i = 0; i < mNDim; ++i) {
        size_t idx_i = 0;
        if (!(mAxes[i].index(pos[i], idx_i))) return false;
        addr += idx_i * mAccu[i];
//         std::cout << "pos[" << i << "] = " << pos[i] << " idx["<<i<<"] = " << idx_i << " ";
    }
//     std::cout << std::endl;
    return true;
}

vector<Axis> HistogramBase::getAxes() const {
    return mAxes;
}

vector<vector<double>> HistogramBase::getTable() const {
    return mPointTable;
}

double HistogramBase::getGridSize() const {
    return mGridSize;
}

size_t HistogramBase::getDimension() const {
    return mNDim;
}

class HistogramValue: public HistogramBase
{
public:
    HistogramValue() {}
    HistogramValue(const vector<Axis>& ax);
    virtual ~HistogramValue() {};
    virtual bool set(const vector<double>& pos, double value = 1.0);
    virtual void fill(double value);
    virtual bool get(const vector<double>& pos, double& value) const;
    friend HistogramValue multiply(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue add(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue minus(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue divide(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue operator+(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue operator-(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue operator*(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue operator*(double x, const HistogramValue& h2);
    friend HistogramValue operator*(const HistogramValue& h1, double x);
    friend HistogramValue operator/(const HistogramValue& h1, const HistogramValue& h2);
    void applyFunction(std::function<double(double)> f);
    vector<double>& getRawData();
    const vector<double>& getRawData() const;
    virtual void writeToFile(const string& filename) const;
    virtual void readFromFile(const string& filename);
    void dump() const;
    void normalize();
    virtual HistogramValue reduceDimension(const vector<size_t> dims) const;
protected:
    vector<double>          mValue;
};

HistogramValue::HistogramValue(const vector<Axis>& ax): HistogramBase(ax) {
    mValue.assign(mGridSize, 0.0);
}

bool HistogramValue::set(const vector<double>& pos, double value) {
    size_t addr = 0;
    bool inGrid = address(pos, addr);
    if (inGrid) {
        mValue[addr] = value;
    }
    return inGrid;
}

void HistogramValue::fill(double value) {
    mValue.assign(mGridSize, value);
}

bool HistogramValue::get(const vector<double>& pos, double& value) const {
    size_t addr = 0;
    bool inGrid = address(pos, addr);
    if (inGrid) {
        value = mValue[addr];
    }
    return inGrid;
}

void HistogramValue::writeToFile(const string& filename) const {
    ofstream ofs_histo(filename.c_str());
    vector<double> pos(mNDim, 0.0);
    double val = 0;
    ofs_histo << "# " << mNDim << '\n';
    for (size_t j = 0; j < mNDim; ++j) {
        ofs_histo << mAxes[j].infoHeader() << '\n';
    }
    ofs_histo.setf(std::ios::fixed);
    ofs_histo << std::setprecision(7);
    for (size_t i = 0; i < mGridSize; ++i) {
        for (size_t j = 0; j < mNDim; ++j) {
            pos[j] = mPointTable[j][i];
            ofs_histo << pos[j] << ' ';
        }
        get(pos, val);
        ofs_histo << val << ' ';
        ofs_histo << '\n';
    }
}

void HistogramValue::readFromFile(const string& filename) {
    ifstream ifs_histo(filename.c_str());
    string line;
    string token{" "};
    vector<string> fields;
    // Parse first line
    std::getline(ifs_histo, line);
    splitString(line, token, fields);
    if (fields[0].compare("#") != 0) {
        std::cerr << "Histogram file reads error!" << std::endl;
        std::abort();
    } else {
        mNDim = std::stoul(fields[1]);
    }
    fields.clear();
    // Parse axes
    mAxes.clear();
    mAxes.resize(mNDim);
    for (size_t i = 0; i < mNDim; ++i) {
        std::getline(ifs_histo, line);
        splitString(line, token, fields);
        double lower, width, upper;
        size_t bins;
        bool periodic = false;
        if (fields[0].compare("#") != 0) {
            std::cerr << "Histogram file reads error!" << std::endl;
            std::abort();
        } else {
            lower = std::stod(fields[1]);
            width = std::stod(fields[2]);
            bins = std::stoul(fields[3]);
            int p = std::stoi(fields[4]);
            upper = lower + double(bins) * width;
            periodic = (p != 0) ? true : false;
            mAxes[i] = Axis(lower, upper, bins, periodic);
        }
        fields.clear();
    }
    // Initialize mAccu
    mAccu.resize(mNDim);
    mAccu[0] = 1;
    mGridSize = 1;
    for (size_t i = 0; i < mNDim; ++i) {
        mAccu[i] = (i == 0) ? 1 : (mAccu[i - 1] * mAxes[i - 1].bin());
        mGridSize *= mAxes[i].bin();
    }
    mValue.assign(mGridSize, 0.0);
    // Initialize table
    fillTable();
    vector<double> pos(mNDim, 0);
    size_t data_count = 0;
    while(std::getline(ifs_histo, line)) {
        splitString(line, token, fields);
        if (fields.empty()) {
            continue;
        }
        if (fields[0].compare("#") != 0) {
            if (fields.size() != (mNDim + 1)) {
                std::cerr << "Histogram file reads error!" << std::endl;
                std::abort();
            }
            double value;
            for (size_t j = 0; j < mNDim; ++j) {
                pos[j] =std::stod(fields[j]);
            }
            value = std::stod(fields[mNDim]);
            set(pos, value);
            ++data_count;
            fields.clear();
        }
    }
    if (data_count != mGridSize) {
        std::cerr << "Histogram file reads error!" << std::endl;
    }
}

void HistogramValue::dump() const {
    using std::cout;
    std::ios_base::fmtflags f(cout.flags());
    vector<double> pos(mNDim, 0.0);
    double val = 0;
    cout.setf(std::ios::fixed);
    cout << std::setprecision(7);
    for (size_t i = 0; i < mGridSize; ++i) {
        for (size_t j = 0; j < mNDim; ++j) {
            pos[j] = mPointTable[j][i];
            cout << pos[j] << ' ';
        }
        get(pos, val);
        for (size_t j = 0; j < mNDim; ++j) {
            cout << val << ' ';
        }
        cout << '\n';
    }
    cout << std::flush;
    cout.flags(f);
}

HistogramValue minus(const HistogramValue& h1, const HistogramValue& h2) {
    HistogramValue h3(h1.mAxes);
    // Assume h1 and h2 have the same axes.
    for (size_t i = 0; i < h3.mGridSize; ++i) {
        h3.mValue[i] = h1.mValue[i] - h2.mValue[i];
    }
    return h3;
}

HistogramValue multiply(const HistogramValue& h1, const HistogramValue& h2) {
    HistogramValue h3(h1.mAxes);
    // Assume h1 and h2 have the same axes.
    for (size_t i = 0; i < h3.mGridSize; ++i) {
        h3.mValue[i] = h1.mValue[i] * h2.mValue[i];
    }
    return h3;
}

HistogramValue add(const HistogramValue& h1, const HistogramValue& h2) {
    HistogramValue h3(h1.mAxes);
    // Assume h1 and h2 have the same axes.
    for (size_t i = 0; i < h3.mGridSize; ++i) {
        h3.mValue[i] = h1.mValue[i] + h2.mValue[i];
    }
    return h3;
}

HistogramValue divide(const HistogramValue& h1, const HistogramValue& h2) {
    HistogramValue h3(h1.mAxes);
    // Assume h1 and h2 have the same axes.
    for (size_t i = 0; i < h3.mGridSize; ++i) {
        if (h2.mValue[i] != 0) {
            h3.mValue[i] = h1.mValue[i] / h2.mValue[i];
        }
    }
    return h3;
}

HistogramValue operator+(const HistogramValue& h1, const HistogramValue& h2) {
    return add(h1, h2);
}

HistogramValue operator-(const HistogramValue& h1, const HistogramValue& h2) {
    return minus(h1, h2);
}

HistogramValue operator*(const HistogramValue& h1, const HistogramValue& h2) {
    return multiply(h1, h2);
}

HistogramValue operator*(double x, const HistogramValue& h2) {
    HistogramValue h3(h2.mAxes);
    for (size_t i = 0; i < h3.mGridSize; ++i) {
        h3.mValue[i] = x * h2.mValue[i];
    }
    return h3;
}

HistogramValue operator*(const HistogramValue& h1, double x) {
    return x * h1;
}

HistogramValue operator/(const HistogramValue& h1, const HistogramValue& h2) {
    return divide(h1, h2);
}

void HistogramValue::applyFunction(std::function<double(double)> f) {
    for (auto it = mValue.begin(); it != mValue.end(); ++it) {
        (*it) = f(*it);
    }
}

vector<double>& HistogramValue::getRawData() {
    return mValue;
}

const vector<double>& HistogramValue::getRawData() const {
    const vector<double>& ret = mValue;
    return ret;
}

void HistogramValue::normalize() {
    double factor = std::accumulate(mValue.begin(), mValue.end(), 0.0);
    if (factor > 0) {
        applyFunction([factor](double x){return x / factor;});
    }
}

HistogramValue HistogramValue::reduceDimension(const vector<size_t> new_dims) const {
    vector<Axis> new_ax;
    for (size_t i = 0; i < new_dims.size(); ++i) {
        new_ax.push_back(mAxes.at(new_dims[i]));
    }
    HistogramValue new_hist(new_ax);
    vector<double> pos(mNDim, 0.0);
    vector<double> new_pos(new_hist.getDimension(), 0.0);
    for (size_t i = 0; i < mGridSize; ++i) {
        double val = 0;
        double new_val = 0;
        for (size_t j = 0; j < mNDim; ++j) {
            pos[j] = mPointTable[j][i];
        }
        for (size_t k = 0; k < new_hist.getDimension(); ++k) {
            new_pos[k] = pos[new_dims[k]];
        }
        get(pos, val);
        new_hist.get(new_pos, new_val);
        new_hist.set(new_pos, new_val + val);
    }
    return new_hist;
}

class ReweightHistogram: public HistogramValue
{
public:
    ReweightHistogram() {}
    ReweightHistogram(const vector<Axis>& ax);
    virtual ~ReweightHistogram() {}
    virtual bool store(const vector<double>& pos, double value = 1.0, double weight = 1.0);
    virtual bool get(const vector<double>& pos, double& value) const;
    size_t getTotalCount() const;
protected:
    vector<double>          mWeightSum;
    vector<size_t>          mCount;
};

ReweightHistogram::ReweightHistogram(const vector<Axis>& ax): HistogramValue(ax) {
    mWeightSum.assign(mGridSize, 0.0);
    mCount.assign(mGridSize, 0);
}

bool ReweightHistogram::get(const vector<double>& pos, double& value) const {
    size_t addr = 0;
    bool inGrid = address(pos, addr);
    if (inGrid) {
        value = mValue[addr];
    }
    return inGrid;
}

bool ReweightHistogram::store(const vector<double>& pos, double value, double weight) {
    size_t addr = 0;
    bool inGrid = address(pos, addr);
    if (inGrid) {
        mValue[addr] += (value * std::exp(weight));
        mWeightSum[addr] += std::exp(weight);
        ++(mCount[addr]);
    }
    return inGrid;
}

size_t ReweightHistogram::getTotalCount() const {
    return std::accumulate(mCount.begin(), mCount.end(), 0);
}

HistogramValue convertToFreeEnergy(const HistogramValue& Punbias, double kbt) {
    HistogramValue FES = Punbias;
    vector<double>& fdata = FES.getRawData();
    bool first_non_zero_value = true;
    double maxValue;
    for (auto& i : fdata) {
        if (i > 0) {
            i = -kbt * std::log(i);
            if (first_non_zero_value) {
                maxValue = i;
                first_non_zero_value = false;
            }
            maxValue = std::max(maxValue, i);
        }
    }
    const vector<double>& pdata = Punbias.getRawData();
    for (size_t i = 0; i < pdata.size(); ++i) {
        if (pdata[i] == 0) {
            fdata[i] = maxValue;
        }
    }
    const double minValue = *std::min_element(fdata.begin(), fdata.end());
    for (auto& i : fdata) {
        i = i - minValue;
    }
    return FES;
}

class HistogramView
{
public:
    HistogramView() {}
    HistogramView(const vector<HistogramValue>& histos);
    HistogramView(const vector<ReweightHistogram>& histos);
    bool addHistogram(HistogramValue& h);
    bool get(const vector<double>& pos, vector<double>& value) const;
    void writeToFile(const string& filename) const;
private:
    vector<reference_wrapper<const HistogramValue>> histogram_references;
    void checkAxes() const;
};

HistogramView::HistogramView (const vector<HistogramValue>& histos):
histogram_references(histos.begin(), histos.end()) {
    checkAxes();
}

HistogramView::HistogramView (const vector<ReweightHistogram>& histos):
histogram_references(histos.begin(), histos.end()) {
    checkAxes();
}

void HistogramView::checkAxes() const {
    if (histogram_references.size() > 0) {
        const vector<Axis> first_ax = histogram_references[0].get().getAxes();
        for (size_t i = 0; i < histogram_references.size(); ++i) {
            const vector<Axis> ax = histogram_references[i].get().getAxes();
            if (ax != first_ax) {
                std::cerr << "Error: Using histograms with different axes in HistogramView class!" << std::endl;
                std::abort();
            }
        }
    }
}

bool HistogramView::addHistogram(HistogramValue& h) {
    histogram_references.push_back(h);
    checkAxes();
    return true;
}

bool HistogramView::get(const vector<double>& pos, vector<double>& value) const {
    value.resize(histogram_references.size(), 0);
    for (size_t i = 0; i < histogram_references.size(); ++i) {
        bool getter_success = histogram_references[i].get().get(pos, value[i]);
        if (!getter_success) {
            return getter_success;
        }
    }
    return true;
}

void HistogramView::writeToFile(const string& filename) const {
    if (histogram_references.size() > 0) {
        ofstream ofs_histo(filename.c_str());
        size_t n_dimension = histogram_references[0].get().getDimension();
        vector<double> pos(n_dimension, 0.0);
        vector<double> val(histogram_references.size(), 0);
        ofs_histo.setf(std::ios::fixed);
        ofs_histo << std::setprecision(7);
        double grid_size = histogram_references[0].get().getGridSize();
        const vector<vector<double>> point_table = histogram_references[0].get().getTable();
        for (size_t i = 0; i < grid_size; ++i) {
            for (size_t j = 0; j < n_dimension; ++j) {
                pos[j] = point_table[j][i];
                ofs_histo << pos[j] << ' ';
            }
            get(pos, val);
            for (size_t j = 0; j < histogram_references.size(); ++j) {
                ofs_histo << val[j] << ' ';
            }
            ofs_histo << '\n';
        }
    }
}

class HistogramFiles : public HistogramBase
{
public:
    HistogramFiles() {}
    HistogramFiles(const vector<Axis>& ax, const string& prefix);
    virtual ~HistogramFiles() {};
    virtual bool store(const vector<double>& pos, const string& data);
    void saveInfo(const string& filename) const;
private:
    string position_to_string(const vector<double>& pos) const;
    string mOutputPrefix;
    vector<ofstream> mFiles;
    vector<string> mFilenames;
};

HistogramFiles::HistogramFiles(const vector<Axis>& ax, const string& prefix): HistogramBase(ax), mOutputPrefix(prefix), mFiles(mGridSize), mFilenames(mGridSize) {
    vector<double> pos(mNDim, 0.0);
    size_t addr = 0;
    for (size_t i = 0; i < mGridSize; ++i) {
        for (size_t j = 0; j < mNDim; ++j) {
            pos[j] = mPointTable[j][i];
        }
        address(pos, addr);
        mFilenames[addr] = mOutputPrefix + position_to_string(pos);
    }
}

string HistogramFiles::position_to_string(const vector<double>& pos) const {
    string result = "_";
    for (size_t i = 0; i < pos.size(); ++i) {
        result += std::to_string(pos[i]) + "_";
    }
    return result;
}

bool HistogramFiles::store(const vector<double>& pos, const string& data) {
    size_t addr = 0;
    bool inGrid = address(pos, addr);
    if (inGrid) {
        mFiles[addr].open(mFilenames[addr].c_str(), ofstream::out | ofstream::app);
        mFiles[addr] << data << '\n';
        mFiles[addr].close();
    }
    return inGrid;
}

void HistogramFiles::saveInfo(const string& filename) const {
    ofstream ofs_info(filename.c_str());
    vector<double> pos(mNDim, 0.0);
    size_t addr = 0;
    for (size_t i = 0; i < mGridSize; ++i) {
        for (size_t j = 0; j < mNDim; ++j) {
            pos[j] = mPointTable[j][i];
        }
        address(pos, addr);
        ofs_info << mFilenames[addr] << " ";
        for (size_t j = 0; j < mNDim; ++j) {
            ofs_info << pos[j] << " ";
        }
        ofs_info << '\n';
    }
}

#endif
