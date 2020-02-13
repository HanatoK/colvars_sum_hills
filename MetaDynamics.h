#ifndef METADYNAMICS_H
#define METADYNAMICS_H

#include "Grid.h"
#include "Axis.h"

#include <vector>
#include <string>
#include <functional>
#include <cmath>
#include <iostream>

#include <json/json.h>

using std::vector;
using std::string;
using std::function;
using std::cout;
using std::ostream;

class Hill
{
public:
    Hill() {};
    Hill(size_t ndims);
    Hill(const vector<double>& centers, const vector<double>& sigmas, double height);
    virtual ~Hill() {}
    void init(size_t ndims);
    virtual void setParameters(const vector<double>& centers, const vector<double>& sigmas, double height);
    virtual double hillEnergy(const vector<double>& pos, const vector<Axis>& axes) const;
    virtual vector<double> hillGradients(const vector<double>& pos, const vector<Axis>& axes) const;
    void debugOutput(ostream& os) const;
    vector<double> mCenters;
    vector<double> mSigmas;
    double mHeight;
};

Hill::Hill(size_t ndims) {
    init(ndims);
}

void Hill::init(size_t ndims) {
    mCenters.resize(ndims);
    mSigmas.resize(ndims);
}

Hill::Hill(const vector<double>& centers, const vector<double>& sigmas, double height) {
    setParameters(centers, sigmas, height);
}

void Hill::setParameters(const vector<double>& centers, const vector<double>& sigmas, double height) {
    mCenters = centers;
    mSigmas = sigmas;
    mHeight = height;
}

double Hill::hillEnergy(const vector<double>& pos, const vector<Axis>& axes) const {
    double result = 0;
    for (size_t i_cv = 0; i_cv < mCenters.size(); ++i_cv) {
        const double dist = axes[i_cv].distance(pos[i_cv], mCenters[i_cv]);
        result += dist * dist / (2.0 * mSigmas[i_cv] * mSigmas[i_cv]);
    }
    result = mHeight * std::exp(-result);
    return result;
}

vector<double> Hill::hillGradients(const vector<double>& pos, const vector<Axis>& axes) const {
    vector<double> grads(pos.size());
    for (size_t i_cv = 0; i_cv < mCenters.size(); ++i_cv) {
        const double dist = axes[i_cv].distance(pos[i_cv], mCenters[i_cv]);
        const double factor = -1.0 * hillEnergy(pos, axes) / (mSigmas[i_cv] * mSigmas[i_cv]);
        grads[i_cv] = dist * factor;
    }
    return grads;
}

void Hill::debugOutput(ostream& os) const {
    std::ios_base::fmtflags f(os.flags());
    os.setf(std::ios::fixed);
    os << std::setprecision(7);
    for (size_t i_cv = 0; i_cv < mCenters.size(); ++i_cv) {
        os << mCenters[i_cv] << "\t ";
    }
    for (size_t i_cv = 0; i_cv < mSigmas.size(); ++i_cv) {
        os << mSigmas[i_cv] << "\t ";
    }
    os << mHeight;
    os << "\n";
    os.flags(f);
}

class MetaDynamics : public HistogramValue
{
public:
    MetaDynamics() {}
    MetaDynamics(const vector<Axis>& ax);
    void addHill(const Hill& hill);
    void compute();
    void readHillsTrajectory(const string& hill_filename);
    void readHillsTrajectory(const Json::Value& obj);
    void readFromJson(const string& json_filename);
    void writeGradients(const string& output_filename) const;
    void writeFakeCount(const string& output_filename) const;
private:
    vector<Hill> mHills;
};

MetaDynamics::MetaDynamics(const vector<Axis>& ax) : HistogramValue(ax), mHills(0) {}

void MetaDynamics::addHill(const Hill& hill) {
    mHills.push_back(hill);
}

void MetaDynamics::compute() {
    for (size_t i_hill = 0; i_hill < mHills.size(); ++i_hill) {
        // TODO: Check if this hill is inside the grid?
        const Hill& h = mHills[i_hill];
        vector<double> pos(mNDim, 0.0);
        for (size_t i_pos = 0; i_pos < mGridSize; ++i_pos) {
            for (size_t j_dim = 0; j_dim < mNDim; ++j_dim) {
                pos[j_dim] = mPointTable[j_dim][i_pos];
            }
            double old_energy = 0;
            get(pos, old_energy);
            double add_energy = h.hillEnergy(pos, mAxes);
            set(pos, old_energy + add_energy);
        }
    }
    this->applyFunction([](double x){return -x;});
}

void MetaDynamics::readHillsTrajectory(const string& hill_filename) {
    ifstream ifs_hilltraj(hill_filename.c_str());
    string line;
    string token{" "};
    vector<string> fields;
    // Parse first line
    std::getline(ifs_hilltraj, line);
    splitString(line, token, fields);
    if (fields[0].compare("#") != 0) {
        std::cerr << "Hills trajectory file reads error!" << std::endl;
        std::abort();
    } else {
        mNDim = std::stoul(fields[1]);
    }
    fields.clear();
    // Parse axes
    mAxes.clear();
    mAxes.resize(mNDim);
    for (size_t i = 0; i < mNDim; ++i) {
        std::getline(ifs_hilltraj, line);
        splitString(line, token, fields);
        double lower, width, upper;
        size_t bins;
        bool periodic = false;
        if (fields[0].compare("#") != 0) {
            std::cerr << "Hills trajectory file reads error!" << std::endl;
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
    while(std::getline(ifs_hilltraj, line)) {
        splitString(line, token, fields);
        if (fields.empty()) {
            continue;
        }
        if (fields.size() < (2*mNDim + 2)) {
            std::cerr << "Hills trajectory file reads error!" << std::endl;
            std::abort();
        }
        if (fields[0].compare("#") != 0) {
            Hill h(mNDim);
            for (size_t i_cv = 0; i_cv < mNDim; ++i_cv) {
                h.mCenters[i_cv] = std::stod(fields[1+i_cv]);
                h.mSigmas[i_cv] = std::stod(fields[1+1*mNDim+i_cv])/2.0;
            }
            h.mHeight = std::stod(fields[1+2*mNDim]);
            addHill(h);
//             cout << line << "\n";
//             h.debugOutput(cout);
        }
        fields.clear();
    }
}

void MetaDynamics::readHillsTrajectory(const Json::Value& obj) {
    // Grid settings
    const Json::Value& grid_config = obj["grid"];
    const unsigned dimension = grid_config.size();
    vector<Axis> a1(dimension);
    for (unsigned i = 0; i < dimension; ++i) {
        double lower = grid_config[i]["lower"].asDouble();
        double upper = grid_config[i]["upper"].asDouble();
        unsigned bins = grid_config[i]["bins"].asUInt();
        bool periodic = false;
        if (grid_config[i].isMember("periodic")) {
            periodic = grid_config[i]["periodic"].asBool();
        }
        a1[i] = Axis(lower, upper, bins, periodic);
    }
    (*this) = MetaDynamics(a1);
    // Output prefix
    const string output_prefix = obj["output"].asString();
    // Trajectories
    const Json::Value& trajectories = obj["trajectory"];
    for (unsigned i = 0; i < trajectories.size(); ++i) {
        const string hill_filename = trajectories[i]["filename"].asString();
        ifstream ifs_hilltraj(hill_filename.c_str());
        string line;
        string token{" "};
        vector<string> fields;
        vector<double> pos(mNDim, 0);
        while(std::getline(ifs_hilltraj, line)) {
            splitString(line, token, fields);
            if (fields.empty() || fields[0].compare("#") == 0) {
                fields.clear();
                continue;
            }
            if (fields.size() < (2*mNDim + 2)) {
                std::cerr << "Hills trajectory file reads error!" << std::endl;
                std::abort();
            }
            if (fields[0].compare("#") != 0) {
                Hill h(mNDim);
                for (size_t i_cv = 0; i_cv < mNDim; ++i_cv) {
                    h.mCenters[i_cv] = std::stod(fields[1+i_cv]);
                    h.mSigmas[i_cv] = std::stod(fields[1+1*mNDim+i_cv])/2.0;
                }
                h.mHeight = std::stod(fields[1+2*mNDim]);
                addHill(h);
            }
            fields.clear();
        }
    }
}

void MetaDynamics::readFromJson(const string& json_filename) {
    ifstream ifs_config(json_filename.c_str());
    Json::CharReaderBuilder reader;
    Json::Value hill_config;
    string json_error;
    Json::parseFromStream(reader, ifs_config, &hill_config, &json_error);
    this->readHillsTrajectory(hill_config);
}

void MetaDynamics::writeGradients(const string& output_filename) const {
    ofstream ofs_grad(output_filename.c_str());
    vector<vector<double>> gradients_table(mNDim, vector<double>(mGridSize, 0.0));
    // compute gradients
    for (size_t i_hill = 0; i_hill < mHills.size(); ++i_hill) {
        // TODO: Check whether hill is inside the grid?
        const Hill& h = mHills[i_hill];
        vector<double> pos(mNDim, 0.0);
        vector<double> grad(mNDim, 0.0);
        for (size_t i_pos = 0; i_pos < mGridSize; ++i_pos) {
            for (size_t j_dim = 0; j_dim < mNDim; ++j_dim) {
                pos[j_dim] = mPointTable[j_dim][i_pos];
            }
            grad = h.hillGradients(pos, mAxes);
            for (size_t j_dim = 0; j_dim < mNDim; ++j_dim) {
                gradients_table[j_dim][i_pos] += grad[j_dim];
            }
        }
    }
    // output
    vector<double> pos(mNDim, 0.0);
    vector<double> grad(mNDim, 0.0);
    ofs_grad.setf(std::ios::fixed);
    ofs_grad << std::setprecision(7);
    ofs_grad << "# " << mNDim << '\n';
    for (size_t j = 0; j < mNDim; ++j) {
        ofs_grad << mAxes[j].infoHeader() << '\n';
    }
    for (size_t i_pos = 0; i_pos < mGridSize; ++i_pos) {
        for (size_t j_dim = 0; j_dim < mNDim; ++j_dim) {
            pos[j_dim] = mPointTable[j_dim][i_pos];
            ofs_grad << pos[j_dim] << ' ';
        }
        for (size_t j_dim = 0; j_dim < mNDim; ++j_dim) {
            grad[j_dim] = gradients_table[j_dim][i_pos];
            ofs_grad << grad[j_dim] << ' ';
        }
        ofs_grad << '\n';
    }
}

void MetaDynamics::writeFakeCount(const string& output_filename) const {
    ofstream ofs_count(output_filename.c_str());
    vector<double> pos(mNDim, 0.0);
    ofs_count.setf(std::ios::fixed);
    ofs_count << std::setprecision(7);
    ofs_count << "# " << mNDim << '\n';
    for (size_t j = 0; j < mNDim; ++j) {
        ofs_count << mAxes[j].infoHeader() << '\n';
    }
    for (size_t i_pos = 0; i_pos < mGridSize; ++i_pos) {
        for (size_t j_dim = 0; j_dim < mNDim; ++j_dim) {
            pos[j_dim] = mPointTable[j_dim][i_pos];
            ofs_count << pos[j_dim] << ' ';
        }
        ofs_count << 1 << '\n';
    }
}

#endif
