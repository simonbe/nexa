#include "NetworkFoldiak.h"

ConnectionModifierFoldiak::ConnectionModifierFoldiak(float eta1, float eta2, float eta3, float alpha, float beta, bool lateral)
{
	m_eta1 = eta1; // f: beta
	m_eta2 = eta2; // f: alpha
	m_eta3 = eta3; // f: -gamma
	m_alpha = alpha; // f: p
	m_lateral = lateral;
	m_beta = beta; // f: lambda

	m_eventId = 8;

	m_transferFunction = new TransferFoldiak(m_beta,m_alpha,m_eta3);
}

void ConnectionModifierFoldiak::Modify()
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
	float dw, ds;

	vector<float> preValues = m_connectionFixed->GetPreValues(m_idsPost[0]);

	for(int j=0;j<postValues.size();j++)
	{	
		if(postValues[j] == 1) // temporary test
		{
			postId = m_idsPost[j];
			vector<long> preIds = m_connectionFixed->GetPreIds(m_idsPost[j]);

			for(int i=0;i<preValues.size();i++)
			{
				//if(preValues[i] == 1) // temporary test
				//{
					preId = preIds[i];//(*preIds)[j][i];

					weight = network()->GetWeight(preId,postId);
					// feedforward weights
					if(m_lateral == false)
					{
						// orig foldiak: m_eta1*postValues[j]*(preValues[i] - weight);
						// Falconbridge modification (also continuous output in transf fcn): m_eta1*postValues[j]*(preValues[i] - postValues[j]*weight);
						dw = m_eta1*postValues[j]*(preValues[i] - weight); // beta

						weight += dw;
					}
					else
					{
						// lateral weights
						if(preId == postId)
							weight=0;
						else
						{
							dw = -m_eta2*(postValues[j]*preValues[i] - m_alpha*m_alpha); // alpha
							weight += dw;

							if(weight>0)
								weight = 0;

							if(weight<-1.0)
								bool b = false;
						}
					}

					network()->SetWeight(weight,preId,postId);
				//}
			}
		}

		// sensitivity
/*		if(m_lateral == false)
		{
			ds = m_eta3*(m_alpha-postValues[j]);
			m_s[j] += ds;

			((RateUnit*)m_connectionFixed->PostLayer()->network()->GetUnitFromId(m_idsPost[j]))->AddInhibBeta(m_s[j]);
		}
*/
	}
}

void ConnectionModifierFoldiak::Simulate(UnitModifier* e)
{
}

void ConnectionModifierFoldiak::SetConnection(Connection* c)
{
	m_connectionFixed = (ConnectionFixed*)c;
}


void ConnectionModifierFoldiak::Initialize(Connection* connection)
{
	network(connection->network());

	vector<float> postValues = m_connectionFixed->GetPostValues();
	m_s = vector<float>(postValues.size());

	m_connectionFixed = connection;
	m_firstRun = true;
	m_idsPost = m_connectionFixed->GetPostIds();
}