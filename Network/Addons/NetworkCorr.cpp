#include "NetworkCorr.h"

void ProjectionModifierPearson::Initialize(Projection* Projection)
{
	m_projection = Projection;

	vector<long> localIds = m_projection->GetPostLocalIds();
	unsigned long postNrUnits = localIds.size();
	unsigned long preNrUnits = m_projection->PreLayer()->GetNrUnitsTotal();

	meanI = vector<float>(preNrUnits,0.0);
	meanJ = vector<float>(postNrUnits,0.0);
	XmeanI = vector<float>(preNrUnits,0.0);
	XmeanJ = vector<float>(postNrUnits,0.0);
	M2I = vector<float>(preNrUnits,0.0);
	M2J = vector<float>(postNrUnits,0.0);
	varianceI = vector<float>(preNrUnits,0.0);
	varianceJ = vector<float>(postNrUnits,0.0);
	variance_nI = vector<float>(preNrUnits,0.0);
	variance_nJ = vector<float>(postNrUnits,0.0);

	rIJ = vector<vector<float> >(preNrUnits,vector<float>(postNrUnits,0.0));
	meanIJ = vector<vector<float> >(preNrUnits,vector<float>(postNrUnits,0.0));

	m_n = 0;
}

void ProjectionModifierPearson::Modify()
{
	float prntaupdt = 0.01;

	vector<float> postValues = m_projection->GetPostValues();
	vector<float> preValues = m_projection->GetPreValues((m_projection->GetPostIds())[0]); // assumes used on a full Projection, should  otherwise loop over all posts

	// online, incremental versions of mean and variance taken from http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance / Knuth (1998), Seminumerical algorithms
	m_n++;

	float delta,deltaX, x,y, corr;
	for(int i=0;i<preValues.size();i++)
	{
		//meanI[i] += (preValues[i] - meanI[i])*prntaupdt;
		x = preValues[i];
		delta = x - meanI[i];
		meanI[i] = meanI[i] + delta/(float)m_n;
		M2I[i] = M2I[i] + delta*(x - meanI[i]); // uses new val of mean

		variance_nI[i] = M2I[i] / (float)m_n; // no need to do every pass
		varianceI[i] = M2I[i] / (float)(m_n-1); // Bessel's correction
	}
	
	for(int j=0;j<postValues.size();j++)
	{
		//meanJ[j] += (postValues[j] - meanJ[j])*prntaupdt;
		x = postValues[j];
		delta = x - meanJ[j];
		meanJ[j] = meanJ[j] + delta/(float)m_n;
		M2J[j] = M2J[j] + delta*(x - meanJ[j]); // uses new val of mean

		variance_nJ[j] = M2J[j] / (float)m_n; // no need to do every pass
		varianceJ[j] = M2J[j] / (float)(m_n-1);

		for(int i=0;i<preValues.size();i++) 
		{
			// own implementation of final expression

			x = (postValues[j] - meanJ[j])*(preValues[i] - meanI[i]);
			deltaX = x - meanIJ[i][j];
			meanIJ[i][j] = meanIJ[i][j] + deltaX/(float)m_n;

			if(fabs(variance_nI[i]) < EPS || fabs(variance_nJ[j]) < EPS)
				corr = 0;
			else
				corr = meanIJ[i][j]/(sqrt(variance_nI[i])*sqrt(variance_nJ[j]));
			
			//y = preValues[i];
			//meanIJ[i][j] += ((y-meanI[i])*(x-meanJ[j]) - meanIJ[i][j]) * prntaupdt;
			//rIJ[i][j] = meanIJ[i][j] / (sqrt(varianceI[i])*sqrt(varianceJ[j]));//( (y-meanI[i])*(x-meanJ[j]) / (sqrt(varianceI[i])*sqrt(varianceJ[j])) - rIJ[i][j] ) * prntaupdt;


			// own definition of corr or zero vectors (group them), and do not group zero with something
			// 1.0 - absolute(pearson corr)
			if(fabs(meanJ[j]) < EPS && fabs(meanI[i]) < EPS)
				rIJ[i][j] = 0.0;
			else if(fabs(meanJ[j]) < EPS || fabs(meanI[i]) < EPS)
				rIJ[i][j] = 1.0;
			else
				rIJ[i][j] = 1.0 - fabs(corr);
		}
	}
}