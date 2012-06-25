#include "NetworkTriesch.h"

ConnectionModifierTriesch::ConnectionModifierTriesch()
{
	m_eventId = 10;
	m_etaHebb =0.05;
	m_beta = 0.2; // go higher over time (eg. -> 0.5, step size: 0.05)

	m_transferFunction = new TransferTriesch(0.005, 1.0, false);
}

ConnectionModifierTriesch::ConnectionModifierTriesch(float etaHebb, float beta, float etaIP, float mu, bool thresholded)
{
	m_eventId = 10;
	m_etaHebb = etaHebb;//0.05;
	m_beta = beta;//0.2; // go higher over time (eg. -> 0.5, step size: 0.05)

	m_transferFunction = new TransferLinear(thresholded,mu,0.001,0.005);//new TransferTriesch(etaIP, mu, thresholded);
}

void ConnectionModifierTriesch::Initialize(Connection* connection)
{
	network(connection->network());

	vector<float> postValues = m_connectionFixed->GetPostValues();
	m_connectionFixed = connection;
	m_firstRun = true;
	m_idsPost = m_connectionFixed->GetPostIds();
}

void ConnectionModifierTriesch::SetConnection(Connection* c)
{
	m_connectionFixed = (ConnectionFixed*)c;
}

void ConnectionModifierTriesch::Modify()
{
	if(IsOn() == false) return;

	int nodeId = m_connectionFixed->PreLayer()->network()->MPIGetNodeId(); // will be put in mpi-class

	if(nodeId == 0) 
	{
		cout<<".";
		cout.flush();
	}

	vector<float> postValues = m_connectionFixed->GetPostValues();
	vector<vector<long> >* preIds;// = m_connectionFixed->PreIds(); // move to initializer (until Invalidate().. called)

	long preId, postId;
	float weight;
	double totWeights;
	float dw, normW;
	float x,y;

	// need max post id for N-function
	vector<long> maxPostIds;// m_connectionFixed->PostLayer()->MPI()->MPIGetMaxIdInHypercolumns(); // usually only one
	bool isMaxValue;

	for(int j=0;j<postValues.size();j++)
	{	
		if(postValues[j] == 1) // temporary test
		{
			vector<float> preValues = m_connectionFixed->GetPreValues(m_idsPost[j]);
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
				totWeights += pow((double)weight,2.0); // fabs(weight); // += weight;
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

float ConnectionModifierTriesch::N(float y, bool isMaxValue)
{
	// standard
	if(isMaxValue)
		return 1.0;
	else 
		return -m_beta;
}

void ConnectionModifierTriesch::Simulate(UnitModifier* e)
{

}

void ConnectionModifierTriesch::SetMu(float mu)
{
	((TransferTriesch*)m_transferFunction)->SetDesiredActivity(mu);
}

void ConnectionModifierTriesch::Clear()
{
	// 1. Clear weight values

	vector<float> postValues = m_connectionFixed->GetPostValues();
	vector<vector<long> >* preIds;// = m_connectionFixed->PreIds(); // move to initializer (until Invalidate().. called)

	long preId, postId;

	// need max post id for N-function
	for(int j=0;j<postValues.size();j++)
	{	
		vector<float> preValues = m_connectionFixed->GetPreValues(m_idsPost[j]);
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