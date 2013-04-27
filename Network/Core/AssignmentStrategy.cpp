#include "AssignmentStrategy.h";
vector<float> SequentialAssignmentStrategy::prepareValues(int idx, const vector<vector<float> >& values) {
	return values[idx];
}

vector<float> ConstLagAssignmentStrategy::prepareValues(int idx, const vector<vector<float> >& values) {
	int interval=values.size()/(_lag+1)-1;
	vector<float> out;
	for(int i=0;i<values[idx].size();i++) {
		for(int j=0;j<=_lag;j++) {
			if(j*interval>idx)
				out.push_back(0);
			else
				out.push_back(values[(idx-interval*j)%values.size()][i]);
		}
	}
	return out;
}

vector<float> ConstLagAssignmentStrategy2::prepareValues(int idx, const vector<vector<float> >& values) {
	int interval=_len/(_lag+1);
	vector<float> out;
	for(int i=0;i<values[idx].size();i++) {
		for(int j=0;j<=_lag;j++) {
			if(j*interval>idx)
				out.push_back(0);
			else
				out.push_back(values[(idx-interval*j)%values.size()][i]);
		}
	}
	return out;
}

vector<float> ExpansionAssignmentStrategy::prepareValues(int idx, const vector<vector<float> >& values) {
	vector<float> out(_n,0.0f);
	for(int i=0;i<values[idx].size();++i)
		out[i]=values[idx][i];
	return out;
}