#include "NetworkTriesch.h"

ProjectionModifierTriesch::ProjectionModifierTriesch()
{
	m_eventId = 10;
	m_etaHebb =0.05;
	m_beta = 0.2; // go higher over time (eg. -> 0.5, step size: 0.05)

	m_transferFunction = new TransferTriesch(0.005, 1.0, false);
}

ProjectionModifierTriesch::ProjectionModifierTriesch(float etaHebb, float beta, float etaIP, float mu, bool thresholded)
{
	m_eventId = 10;
	m_etaHebb = etaHebb;//0.05;
	m_beta = beta;//0.2; // go higher over time (eg. -> 0.5, step size: 0.05)

	m_transferFunction = new TransferLinear(thresholded,mu,0.001,0.005);//new TransferTriesch(etaIP, mu, thresholded);
}

void ProjectionModifierTriesch::Initialize(Projection* Projection)
{
	network(Projection->network());

	vector<float> postValues = m_projectionFixed->GetPostValues();
	m_projectionFixed = Projection;
	m_firstRun = true;
	m_idsPost = m_projectionFixed->GetPostIds();
}

void ProjectionModifierTriesch::SetProjection(Projection* c)
{
	m_projectionFixed = (ProjectionFixed*)c;
}

void ProjectionModifierTriesch::Modify()
{
	if(IsOn() == false) return;

	int processId = m_projectionFixed->PreLayer()->network()->MPIGetNodeId(); // will be put in mpi-class

	if(processId == 0) 
	{
		cout<<".";
		cout.flush();
	}

	vector<float> postValues = m_projectionFixed->GetPostValues();
	vector<vector<long> >* preIds;// = m_projectionFixed->PreIds(); // move to initializer (until Invalidate().. called)

	long preId, postId;
	float weight;
	double totWeights;
	float dw, normW;
	float x,y;

	// need max post id for N-function
	vector<long> maxPostIds;// m_projectionFixed->PostLayer()->MPI()->MPIGetMaxIdInHypercolumns(); // usually only one
	bool isMaxValue;

	for(int j=0;j<postValues.size();j++)
	{	
		if(postValues[j] == 1) // temporary test
		{
			vector<float> preValues = m_projectionFixed->GetPreValues(m_idsPost[j]);
			postId = m_idsPost[j];
			y = postValues[j];

			isMaxValue = false;

			for(int k=0;k<maxPostIds.size();k++)
			{
				if(postId == maxPostIds[k])
				{
					isMaxValue = true;
					break;
				}
			}

			float Ny = N(y,isMaxValue);

			for(int i=0;i<preValues.size();i++)
			{
				preId = (*preIds)[j][i];
				x = preValues[i];

				weight = network()->GetWeight(preId,postId);

				dw = x*y*Ny; // (alt: N implemented as an outside WTA)
				weight = weight + m_etaHebb*dw;

				network()->SetWeight(weight,preId,postId);
			}

			// weight normalization

			totWeights = 0;

			for(int i=0;i<preValues.size();i++)
			{
				preId = (*preIds)[j][i];

				weight = network()->GetWeight(preId,postId);
				totWeights += weight*weight;//pow((double)weight,2.0); // fabs(weight); // += weight;
			}

			totWeights = sqrt(totWeights);

			if(totWeights > 0.0)
				for(int i=0;i<preValues.size();i++)
				{
					preId = (*preIds)[j][i];

					weight = network()->GetWeight(preId,postId);
					weight = weight/totWeights;
					network()->SetWeight(weight,preId,postId);
				}
		}
	}
}

float ProjectionModifierTriesch::N(float y, bool isMaxValue)
{
	// standard
	if(isMaxValue)
		return 1.0;
	else 
		return -m_beta;
}

void ProjectionModifierTriesch::Simulate(UnitModifier* e)
{

}

void ProjectionModifierTriesch::SetMu(float mu)
{
	((TransferTriesch*)m_transferFunction)->SetDesiredActivity(mu);
}

void ProjectionModifierTriesch::Clear()
{
	// 1. Clear weight values

	vector<float> postValues = m_projectionFixed->GetPostValues();
	vector<vector<long> >* preIds;// = m_projectionFixed->PreIds(); // move to initializer (until Invalidate().. called)

	long preId, postId;

	// need max post id for N-function
	for(int j=0;j<postValues.size();j++)
	{	
		vector<float> preValues = m_projectionFixed->GetPreValues(m_idsPost[j]);
		postId = m_idsPost[j];
		for(int i=0;i<preValues.size();i++)
		{
			preId = (*preIds)[j][i];

			network()->SetWeight(0.0,preId,postId);
		}
	}

	// 2. Clear transfer function values
	((TransferTriesch*)m_transferFunction)->Clear();
}