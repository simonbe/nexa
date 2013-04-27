#pragma once
#include "Network.h"
class IAssignmentStrategy {
public:
	virtual vector<float> prepareValues(int, const vector<vector<float> >&) = 0;
};

class SequentialAssignmentStrategy:public IAssignmentStrategy {
public:
	vector<float> prepareValues(int idx, const vector<vector<float> >& values);
};

class ExpansionAssignmentStrategy:public IAssignmentStrategy {
public:
	vector<float> prepareValues(int idx, const vector<vector<float> >& values);
	ExpansionAssignmentStrategy(int n) {_n=n;};
private:
	int _n;
};

class ConstLagAssignmentStrategy:public IAssignmentStrategy {
public:
	vector<float> prepareValues(int idx, const vector<vector<float> >& values);
	ConstLagAssignmentStrategy(int lag) {_lag=lag;};
private:
	int _lag;
};

class ConstLagAssignmentStrategy2:public IAssignmentStrategy {
public:
	vector<float> prepareValues(int idx, const vector<vector<float> >& values);
	ConstLagAssignmentStrategy2(int len, int lag) {_len=len;_lag=lag;};
private:
	int _lag;
	int _len;
};