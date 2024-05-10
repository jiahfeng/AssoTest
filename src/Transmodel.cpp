// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <random>
#include <fstream>
#include <string>
#include <sstream>
//#include <Eigen/Dense>
//#include <Eigen/Core>
//#include <Eigen/LU>

using namespace Rcpp;
using namespace Eigen;
using namespace std;


class SurvData {
public:

  VectorXd Y_, CI_;
  MatrixXd X_;

};

class StepFunc {
public:

  VectorXd time_, jump_, cum_;

  void inputjump(const VectorXd& t, const VectorXd& jmp);
  void calcum();
  double eval(double x);
};

class Newton {
protected:
  double tol_;

  VectorXd param_;

public:

  virtual double likelihood() = 0;
  virtual VectorXd score() = 0;
  virtual MatrixXd hess() = 0;
  void iterate();
  void maximize();
};

class CoxRegress : public Newton {
public:

  int n_, p_;
  StepFunc hazard_;
  SurvData data_;
  VectorXd beta_, frailty_, S0_;
  vector<VectorXd> S1_;
  vector<MatrixXd> S2_;

  void inputdata(const VectorXd& Y, const VectorXd& CI, const MatrixXd& X);
  void initialize();
  void setfrailty(const VectorXd& w);
  virtual double likelihood();
  void calS();
  virtual VectorXd score();
  virtual MatrixXd hess();
  void updatehazard();
  VectorXd getbeta() const;
  StepFunc gethazard();
};

class TransCox : public CoxRegress {
public:
  double rho_;
  double (*dG_)(double, double);
  double (*ddG_)(double, double);
  double (*dddG_)(double, double);

  void settransform(double (*dG)(double, double), double (*ddG)(double, double), double (*dddG)(double, double));
  void EM();
};

void StepFunc::inputjump(const VectorXd& t, const VectorXd& jmp) {
  time_ = t;
  jump_ = jmp;
}

void StepFunc::calcum() {
  size_t size = time_.size();
  VectorXd cum = VectorXd::Constant(size, 0);

  for (size_t i = 0; i < size; ++i) {
    if (i > 0) {
      if (time_(i - 1) == time_(i)) {
        cum(i) = cum(i - 1);
      }
      else {
        double temp = 0;
        int j = i;
        do {
          temp += jump_(j);
          ++j;
          if (j == size) break;
        } while (time_(j) == time_(j - 1));
        cum(i) = cum(i - 1) + temp;
      }
    }
    else {
      double temp = 0;
      int j = i;
      do {
        temp += jump_(j);
        ++j;
        if (j == size) break;
      } while (time_(j) == time_(j - 1));
      cum(i) = temp;
    }
  }
  cum_ = cum;
}

double StepFunc::eval(double x) {
  for (size_t i = 0; i < time_.size() - 1; ++i) {
    if (time_(i) <= x && time_(i + 1) > x) return cum_(i);
    else if (i == time_.size() - 1) return cum_(i);
  }
}

void Newton::iterate() {
  VectorXd step;

  MatrixXd hess_inv = hess().inverse();
  step = hess_inv * score();

  param_ = param_ - step;
}

void Newton::maximize() {
  VectorXd step;
  size_t iter = 0;
  tol_ = 1e-5;

  do {
    ++iter;
    MatrixXd hess_inv = hess().inverse();
    step = hess_inv * score();

    param_ = param_ - step;
  } while (step.norm() > tol_ && iter < 300);
}

void CoxRegress::inputdata(const VectorXd& Y, const VectorXd& CI, const MatrixXd& X) {
  data_.Y_ = Y;
  data_.CI_ = CI;
  data_.X_ = X;
  n_ = data_.Y_.size();
  p_ = data_.X_.cols();
}

void CoxRegress::initialize() {
  param_.resize(p_);
  param_.setZero();
  beta_ = param_;
  frailty_ = VectorXd::Constant(data_.Y_.size(), 1);
  VectorXd jumps(data_.Y_.size());
  for (size_t i = 0; i < data_.Y_.size(); ++i) jumps(i) = data_.CI_(i) / data_.CI_.sum();
  hazard_.inputjump(data_.Y_, jumps);
  hazard_.calcum();
}

void CoxRegress::setfrailty(const VectorXd& w) {
  frailty_ = w;
}

void CoxRegress::calS() {
  VectorXd eta = data_.X_ * beta_;
  S0_.resize(n_);
  S1_.resize(n_);
  S2_.resize(n_);
  S0_(0) = (frailty_.array() * eta.array().exp()).sum();
  S1_[0] = VectorXd::Zero(p_);
  S2_[0] = MatrixXd::Zero(p_, p_);
  for (size_t i = 0; i < n_; ++i) {
    S1_[0] = S1_[0] + frailty_(i) * exp(eta(i)) * data_.X_.row(i).transpose();
    S2_[0] = S2_[0] + frailty_(i) * exp(eta(i)) * data_.X_.row(i).transpose() * data_.X_.row(i);
  }
  for (size_t i = 1; i < n_; ++i) {
    S0_(i) = S0_(i - 1) - frailty_(i - 1) * exp(eta(i - 1));
    S1_[i] = S1_[i - 1] - frailty_(i - 1) * exp(eta(i - 1)) * data_.X_.row(i - 1).transpose();
    S2_[i] = S2_[i - 1] - frailty_(i - 1) * exp(eta(i - 1)) * data_.X_.row(i - 1).transpose() * data_.X_.row(i - 1);
  }
  for (size_t i = 1; i < n_; ++i) {
    if (data_.Y_(i) == data_.Y_(i - 1)) {
      S0_(i) = S0_(i - 1);
      S1_[i] = S1_[i - 1];
      S2_[i] = S2_[i - 1];
    }
  }
}

double CoxRegress::likelihood() {
  beta_ = param_;
  VectorXd eta = data_.X_ * beta_;
  double like = 0;
  for (size_t i = 0; i < n_; ++i) {
    if (data_.CI_(i) == 1) {
      like += eta(i) - log(S0_(i));
    }
  }
  return like;
}

VectorXd CoxRegress::score() {
  beta_ = param_;
  VectorXd eta = data_.X_ * beta_;
  VectorXd sc = VectorXd::Zero(p_);

  for (size_t i = 0; i < n_; ++i) {
    if (data_.CI_(i) == 1) {
      sc += data_.X_.row(i).transpose() - S1_[i] / S0_(i);
    }
  }
  return sc;
}

MatrixXd CoxRegress::hess() {
  beta_ = param_;
  VectorXd eta = data_.X_ * beta_;
  MatrixXd hs = MatrixXd::Zero(p_, p_);

  for (size_t i = 0; i < n_; ++i) {
    if (data_.CI_(i) == 1) {
      hs -= S2_[i] / S0_(i) - S1_[i] * S1_[i].transpose() / pow(S0_(i), 2);
    }
  }
  return hs;
}

void CoxRegress::updatehazard() {
  beta_ = param_;
  VectorXd eta = data_.X_ * beta_;
  S0_.resize(n_);
  S0_(0) = (frailty_.array() * eta.array().exp()).sum();
  for (size_t i = 1; i < n_; ++i) {
    S0_(i) = S0_(i - 1) - frailty_(i - 1) * exp(eta(i - 1));
  }
  for (size_t i = 1; i < n_; ++i) {
    if (data_.Y_(i) == data_.Y_(i - 1)) {
      S0_(i) = S0_(i - 1);
    }
  }
  VectorXd jumps(n_);
  for (size_t i = 0; i < n_; ++i) {
    jumps(i) = data_.CI_(i) / S0_(i);
  }
  hazard_.inputjump(data_.Y_, jumps);
  hazard_.calcum();
}

VectorXd CoxRegress::getbeta() const {
  return param_;
}

StepFunc CoxRegress::gethazard() {
  return hazard_;
}

void TransCox::settransform(double (*dG)(double, double), double (*ddG)(double, double), double (*dddG)(double, double)) {
  dG_ = *dG;
  ddG_ = *ddG;
  dddG_ = *dddG;
}

void TransCox::EM() {
  CoxRegress coxtemp;

  coxtemp.inputdata(data_.Y_, data_.CI_, data_.X_);
  coxtemp.initialize();

  VectorXd beta, eta, betadiff, hazarddiff, xi(data_.Y_.size()); //xi is the frailty
  StepFunc hazard;
  double diff;
  do {
    beta = coxtemp.getbeta();
    eta = data_.X_ * beta;
    hazard = coxtemp.gethazard();

    //E step
    for (size_t i = 0; i < data_.Y_.size(); ++i) {
      if (data_.CI_(i) == 0) {
        xi(i) = (*dG_)(hazard.cum_(i) * exp(eta(i)), rho_);
      }
      else {
        xi(i) = -(*ddG_)(hazard.cum_(i) * exp(eta(i)), rho_) / (*dG_)(hazard.cum_(i) * exp(eta(i)), rho_) + (*dG_)(hazard.cum_(i) * exp(eta(i)), rho_);
      }
    }

    //M step
    coxtemp.setfrailty(xi);
    coxtemp.calS();
    cout << coxtemp.likelihood() << endl;
    coxtemp.iterate();
    coxtemp.updatehazard();

    betadiff = coxtemp.getbeta() - beta;
    betadiff = betadiff.array().abs();
    hazarddiff = coxtemp.gethazard().jump_ - hazard.jump_;
    hazarddiff = hazarddiff.array().abs();
    if (p_ > 0) diff = betadiff.maxCoeff(); else diff = 0;
  } while (diff > 1e-8 || hazarddiff.sum() > 1e-8);

  beta = coxtemp.getbeta();
  eta = data_.X_ * beta;
  param_ = coxtemp.getbeta();
  hazard_ = coxtemp.gethazard();
  for (size_t i = 0; i < data_.Y_.size(); ++i) {
    if (data_.CI_(i) == 0) {
      xi(i) = (*dG_)(hazard_.cum_(i) * exp(eta(i)), rho_);
    }
    else {
      xi(i) = -(*ddG_)(hazard_.cum_(i) * exp(eta(i)), rho_) / (*dG_)(hazard_.cum_(i) * exp(eta(i)), rho_) + (*dG_)(hazard_.cum_(i) * exp(eta(i)), rho_);
    }
  }
  frailty_ = xi;
}

void inputXY(VectorXd& Y, VectorXd& CI, MatrixXd& X) {
  std::string datafile;
  datafile = "data.csv";
  std::ifstream indata;
  std::string line;

  indata.open(datafile);
  int n = 5000, p = 3;

  Y.resize(n);
  CI.resize(n);
  X.setZero(n, p);
  if (indata.is_open()) {
    int obs = 0;
    while (getline(indata, line)) {
      string out;
      std::istringstream stream(line);
      getline(stream, out, ',');
      Y(obs) = atof(out.c_str());
      getline(stream, out, ',');
      CI(obs) = atof(out.c_str());
      for (size_t i = 0; i < p; ++i) {
        getline(stream, out, ',');
        X(obs, i) = atof(out.c_str());
      }
      obs++;
    }
    n = obs;
  }
  indata.close();

  MatrixXd Xtmp = X.topRows(n);
  VectorXd CItmp = CI.head(n);
  VectorXd Ytmp = Y.head(n);
  Y = Ytmp;
  X = Xtmp;
  CI = CItmp;
}

double dG(double x, double rho) {
  if (rho == 1) {
    return 1;
  }
  else {
    return 1 / (1 + x);
  }
}

double ddG(double x, double rho) {
  if (rho == 1) {
    return 0;
  }
  else {
    return -pow(1 + x, -2);
  }
}

double dddG(double x, double rho) {
  if (rho == 1) {
    return 0;
  }
  else {
    return 2 * pow(1 + x, -3);
  }
}

//[[Rcpp::export]]
Rcpp::List Trans_EM_rcpp(const Eigen::VectorXd Y, const Eigen::VectorXd CI, const Eigen::MatrixXd X, double RHO) {
  int p_ = X.cols();

  CoxRegress coxtemp;

  coxtemp.inputdata(Y, CI, X);
  coxtemp.initialize();

  VectorXd beta, eta, betadiff, hazarddiff, xi(Y.size()); //xi is the frailty
  StepFunc hazard;
  double diff;
  do {
    beta = coxtemp.getbeta();
    eta = X * beta;
    hazard = coxtemp.gethazard();

    //E step
    for (size_t i = 0; i < Y.size(); ++i) {
      if (CI(i) == 0) {
        xi(i) = dG(hazard.cum_(i) * exp(eta(i)), RHO);
      }
      else {
        xi(i) = -ddG(hazard.cum_(i) * exp(eta(i)), RHO) / dG(hazard.cum_(i) * exp(eta(i)), RHO) + dG(hazard.cum_(i) * exp(eta(i)), RHO);
      }
    }

    //M step
    coxtemp.setfrailty(xi);
    coxtemp.calS();
    coxtemp.iterate();
    coxtemp.updatehazard();

    betadiff = coxtemp.getbeta() - beta;
    betadiff = betadiff.array().abs();
    hazarddiff = coxtemp.gethazard().jump_ - hazard.jump_;
    hazarddiff = hazarddiff.array().abs();
    if (p_ > 0) diff = betadiff.maxCoeff(); else diff = 0;
  } while (diff > 1e-8 || hazarddiff.sum() > 1e-8);


  return Rcpp::List::create(Rcpp::Named("beta") = coxtemp.getbeta(), Rcpp::Named("hazard") = coxtemp.gethazard().jump_);
}

int main(int argc, char* argv[]) {
  MatrixXd X;
  VectorXd Y, CI;

  inputXY(Y, CI, X);

  TransCox cox;
  cox.inputdata(Y, CI, X);
  cox.initialize();
  cox.settransform(dG, ddG, dddG);
  cox.EM();

  cout << cox.getbeta().transpose() << endl;

  return 0;
}
