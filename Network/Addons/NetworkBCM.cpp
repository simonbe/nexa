#include "NetworkBCM.h"

ConnectionModifierBCM::ConnectionModifierBCM()
{
	m_eventId = 13;
	m_eta = 0.1;
	m_decay = 0.5;
	m_tau = 10;

	m_transferFunction = new TransferSigmoid(false);
}

ConnectionModifierBCM::ConnectionModifierBCM(float eta, float decay, float tau)
{
	m_eventId = 13;
	m_eta = eta;
	m_decay = decay;
	m_tau = tau;

	m_transferFunction = new TransferSigmoid(false);
}

void ConnectionModifierBCM::SetEta(float eta)
{
	m_eta = eta;
}

void ConnectionModifierBCM::Initialize(Connection* connection)
{
	network(connection->network());

	m_connectionFixed = connection;
	m_firstRun = true;
	m_idsPost = m_connectionFixed->GetPostIds();
}

void ConnectionModifierBCM::SetConnection(Connection* c)
{
	m_connectionFixed = (ConnectionFixed*)c;
}

void ConnectionModifierBCM::Modify()
{
	if(IsOn() == false) return;

	int nodeId = m_connectionFixed->PreLayer()->network()->MPIGetNodeId(); // will be put in mpi-class

	if(nodeId == 0) 
	{
		cout<<".";
		cout.flush();
	}

	vector<float> postValues = m_connectionFixed->GetPostValues();

	// this does not work correctly atm
	//vector<vector<long> >* preIds = m_connectionFixed->PreIds(); // move to initializer (until Invalidate().. called)

	if(m_firstRun == true)
	{
		m_thresholds = vector<float>(postValues.size(),0.1);
		m_firstRun = false;
	}

	long preId, postId;
	float weight;
	double totWeights;
	float dw, normW;
	float x,y;
	float phi;

	for(int j=0;j<postValues.size();j++)
	{
		vector<float> preValues = m_connectionFixed->GetPreValues(m_idsPost[j]);
		vector<long> preIds = m_connectionFixed->GetPreIds(m_idsPost[j]);
		postId = m_idsPost[j];
		y = postValues[j];

//		if(m_thresholds[j] == 0.0)
//			phi = 0.0; // 1?
//		else
			phi = m_eta*y*( y - m_thresholds[j] ) / m_thresholds[j]; // law and cooper bcm

		for(int i=0;i<preValues.size();i++)
		{
			preId = preIds[i];//(*preIds)[j][i];
			x = preValues[i];

			if(preId == postId)
				bool b = false;

			weight = network()->GetWeight(preId,postId);
			
			weight = weight + phi*x - m_eta*m_decay*weight;

			network()->SetWeight(weight,preId,postId);
		}

		m_thresholds[j] = m_thresholds[j] + (y*y-m_thresholds[j])/m_tau;

		// weight normalization

		/*totWeights = 0;

		for(int i=0;i<preValues.size();i++)
		{
			preId = preIds[i];//(*preIds)[j][i];

			weight = network()->GetWeight(preId,postId);
			totWeights += pow((double)weight,2.0); // fabs(weight); // += weight;
		}

		totWeights = sqrt(totWeights);

		if(totWeights > 0.0)
		for(int i=0;i<preValues.size();i++)
		{
			preId = preIds[i];//(*preIds)[j][i];

			weight = network()->GetWeight(preId,postId);
			weight = weight/totWeights;
			network()->SetWeight(weight,preId,postId);
		}
		*/
	}
}

void ConnectionModifierBCM::Simulate(UnitModifier* e)
{
}

/*
		elif learning_rule==1: # bcm
            count=0
            for n from 0<=n < num_neurons:
                phi[n]=eta*y[n]*(y[n]-th[n])
                for i from 0<=i<num_inputs:
                    w[count]=w[count]+phi[n]*x[count]-eta*decay*w[count]
                        
                    count=count+1
                    
                th[n]=th[n]+(y[n]*y[n]-th[n])/tau
            
        elif learning_rule==2: # law and cooper bcm
            count=0
            for n from 0<=n < num_neurons:
                phi[n]=eta*y[n]*(y[n]-th[n])/th[n]
                
                for i from 0<=i<num_inputs:
                    w[count]=w[count]+phi[n]*x[count]-eta*decay*w[count]
                    count=count+1
                    
                th[n]=th[n]+(y[n]*y[n]-th[n])/tau
*/