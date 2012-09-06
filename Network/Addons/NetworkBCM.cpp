#include "NetworkBCM.h"

ProjectionModifierBCM::ProjectionModifierBCM()
{
	m_eventId = 13;
	m_eta = 0.1;
	m_decay = 0.5;
	m_tau = 10;

	m_transferFunction = new TransferSigmoid(false);
}

ProjectionModifierBCM::ProjectionModifierBCM(float eta, float decay, float tau)
{
	m_eventId = 13;
	m_eta = eta;
	m_decay = decay;
	m_tau = tau;

	m_transferFunction = new TransferSigmoid(false);
}

void ProjectionModifierBCM::SetEta(float eta)
{
	m_eta = eta;
}

void ProjectionModifierBCM::Initialize(Projection* Projection)
{
	network(Projection->network());

	m_projectionFixed = Projection;
	m_firstRun = true;
	m_idsPost = m_projectionFixed->GetPostIds();
}

void ProjectionModifierBCM::SetProjection(Projection* c)
{
	m_projectionFixed = (ProjectionFixed*)c;
}

void ProjectionModifierBCM::Modify()
{
	if(IsOn() == false) return;

	int processId = m_projectionFixed->PreLayer()->network()->MPIGetNodeId(); // will be put in mpi-class

	if(processId == 0) 
	{
		cout<<".";
		cout.flush();
	}

	vector<float> postValues = m_projectionFixed->GetPostValues();

	// this may not work correctly atm
	//vector<vector<long> >* preIds = m_projectionFixed->PreIds(); // move to initializer (until Invalidate().. called)

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
		vector<float> preValues = m_projectionFixed->GetPreValues(m_idsPost[j]);
		vector<long> preIds = m_projectionFixed->GetPreIds(m_idsPost[j]);
		postId = m_idsPost[j];
		y = postValues[j];

//		if(m_thresholds[j] == 0.0)
//			phi = 0.0; // or 1
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

		// weight normalization not used

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

void ProjectionModifierBCM::Simulate(UnitModifier* e)
{
}