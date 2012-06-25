#include "NetworkCL.h"

ConnectionModifierCL::ConnectionModifierCL(int sizeOutLayer, float eta, float probB, float biasC)
{
	m_eta = eta;//0.01;
	B = probB;// 0.001;
	C = biasC;// 10
	N = (float)sizeOutLayer;

	m_eventId = 6;

	m_transferFunction = new TransferLinear(false);
}

void ConnectionModifierCL::Initialize(Connection* connection)
{
	network(connection->network());

	vector<float> postValues = m_connectionFixed->GetPostValues();
	m_p = vector<float>(postValues.size());
	m_b = vector<float>(postValues.size());

	m_connectionFixed = connection;
	m_firstRun = true;
	m_idsPost = m_connectionFixed->GetPostIds();
}

void ConnectionModifierCL::SetConnection(Connection* c)
{
	m_connectionFixed = (ConnectionFixed*)c;
}


void ConnectionModifierCL::SetC(float value)
{
	C = value;
}

void ConnectionModifierCL::Modify()
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
	float dw;
	
	

	for(int j=0;j<postValues.size();j++)
	{	
		if( postValues[j] == 1.0) // assume WTA(1.0) is used
		{
			postId = m_idsPost[j];

			vector<float> preValues = m_connectionFixed->GetPreValues(m_idsPost[j]);

			for(int i=0;i<preValues.size();i++)
			{
				preId = (*preIds)[j][i];

				weight = network()->GetWeight(preId,postId);
				dw = m_eta*(preValues[i] - weight);

				weight = weight + dw;

				network()->SetWeight(weight,preId,postId);
			}

			// weight normalization

			totWeights = 0;
			
			for(int i=0;i<preValues.size();i++)
			{
				preId = (*preIds)[j][i];

				weight = network()->GetWeight(preId,postId);
				totWeights += weight;//pow((double)weight,2.0);
			}

			//totWeights = sqrt(totWeights);

			for(int i=0;i<preValues.size();i++)
			{
				preId = (*preIds)[j][i];

				weight = network()->GetWeight(preId,postId);
				weight = weight/totWeights;
				network()->SetWeight(weight,preId,postId);
			}
		}
		
		// threshold ("conscience" implementation from DeSieno)

		m_p[j] = m_p[j] + B*(postValues[j] - m_p[j]);
		m_b[j] = C*(1/N - m_p[j]);

		((RateUnit*)m_connectionFixed->PostLayer()->network()->GetUnitFromId(m_idsPost[j]))->AddInhibBeta(m_b[j]);
	}
}

void ConnectionModifierCL::Simulate(UnitModifier* e)
{
}



// Invalidate()...