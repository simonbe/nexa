#include "NetworkSTDP.h"

/*function [mod_weight,exp_deca1,exp_deca2] = stdpm(pre_synap,post_synap,weight,preexp_decay1,preexp_decay2,dw)
    Wmax=10000; %12000
    Wmin=0;

    exp_deca1 = synapse_stdp(preexp_decay1,dw,pre_synap);
    exp_deca2 = synapse_stdp(preexp_decay2,dw,post_synap);

    if (pre_synap == 1)
        
        weight=weight+exp_deca2;%1
%         exp_deca2=preexp_decay2;
    elseif (post_synap == 1)
        
        weight=weight-exp_deca1;%2
%         exp_deca1=preexp_decay1;

%         exp_deca1=preexp_decay1;
%         exp_deca2=preexp_decay2;
    end
    
    if weight>Wmax
        weight=Wmax;
    end
    if weight<Wmin
        weight=Wmin;
    end
     
    mod_weight=weight;
    */


ProjectionModifierSTDP::ProjectionModifierSTDP(bool inhibitoryWeights)
{
	m_eventId = 15;
	//m_transferFunction = new TransferLinear(false);

	if(inhibitoryWeights == true)
		m_scaleFactor = -1;
	else
		m_scaleFactor = 1;

	Wmax=10000;//m_scaleFactor*10000; // 12000
    Wmin=0;//m_scaleFactor*0;
	Wadj=5;
}

void ProjectionModifierSTDP::Initialize(Projection* Projection)
{
	network(Projection->network());

	m_projectionFixed = Projection;
	m_firstRun = true;
	m_idsPost = m_projectionFixed->GetPostIds();
}

void ProjectionModifierSTDP::Modify()
{
	if(!IsOn()) return;
	// Should iterate over GetPostSpikes, GetPreSpikes.. etc
	
	// currently a fix to this
	if(this->m_network->MPIGetNrProcs()>1)
		this->m_projectionFixed->PreLayer()->MPI()->MPIMakeLayerValuesLocal();	

	vector<float> postValues = m_projectionFixed->GetPostValues();
	if(m_prevPostValues.size()==0)
		m_prevPostValues = postValues;

	//vector<vector<long> >* preIds = m_projectionFixed->GetPreIdsAll();// = m_projectionFixed->PreIds(); // move to initializer (until Invalidate().. called)
	vector<long>* preIds = m_projectionFixed->GetPreIdsAll();// = m_projectionFixed->PreIds(); // move to initializer (until Invalidate().. called)

	//vector<long> postIds = m_projectionFixed->GetPostIds();

	if(m_firstRun == true)
	{
		// initialize vectors
		//preexp_decay2 = vector<float>(postValues.size(),0.0);
		preexp_decay2 = vector<vector<float> >(postValues.size());
		preexp_decay1 = vector<vector<float> >(postValues.size());
		for(int j=0;j<postValues.size();j++)
		{
			vector<long> preIds = m_projectionFixed->GetPreIds(m_idsPost[j]);
			//vector<float> preValues = m_projectionFixed->GetPreValues(m_idsPost[j]);
			preexp_decay2[j] = vector<float>(preIds.size(),0.0);//preValues.size(),0.0);
			preexp_decay1[j] = vector<float>(preIds.size(),0.0);
		}

		m_firstRun = false;
	}

	long preId, postId;
	float weight;
	double totWeights;
	float dw, normW;
	float x,y;

	//exp_deca1 = synapse_stdp(preexp_decay1,dw,pre_synap);
	//exp_deca2 = synapse_stdp(preexp_decay2,dw,post_synap);

	for(int j=0;j<postValues.size();j++)
	{
		vector<float> preValues = m_projectionFixed->GetPreValues(m_idsPost[j]);
		postId = m_idsPost[j];

		vector<long> preIds = m_projectionFixed->GetPreIds(m_idsPost[j]);
		
		for(int i=0;i<preValues.size();i++)
		{
			preexp_decay2[j][i] = synapse_stdp(preexp_decay2[j][i],Wadj,m_prevPostValues[j]);//postValues[j]);
			preexp_decay1[j][i] = synapse_stdp(preexp_decay1[j][i],Wadj,preValues[i]);

			bool update = false;
			
			if(preValues[i] > 0.9) // spike
			{
				preId = preIds[i];//(*preIds)[j][i];
				weight = m_network->GetWeight(preId,postId);
				
				weight = weight + m_scaleFactor*preexp_decay2[j][i];
				update = true;
			}
			else if(m_prevPostValues[j] > 0.9)//postValues[j] > 0.9) // spike
			{
				preId = preIds[i];//(*preIds)[j][i];
				weight = m_network->GetWeight(preId,postId);
				float exp_deca1 = 0;
				weight = weight - m_scaleFactor*preexp_decay1[j][i];
				update = true;
			}

			if(update)
			{
				if(m_scaleFactor*weight>Wmax)
					weight=m_scaleFactor*Wmax;
				if(m_scaleFactor*weight<Wmin)
					weight=m_scaleFactor*Wmin;

				m_network->SetWeight(weight,preId,postId);
			}
		}
	}

	m_prevPostValues = postValues;
}


/*
function current = synapse_stdp (I,w,input)

% Ke = 0.0004;
% Ke=0.001;
%Ke=0.4;%0.054;
Ke=0.27;
    if (input == 1)
        current = I-I*Ke+w;
    else if (input == 0)
            current = I-I*Ke;
        end
    end
*/


float ProjectionModifierSTDP::synapse_stdp(float I, float w, float input)
{
/*
% input: spikes
% w: weight
% I: previous potential
% current: updated potential
*/
	float current;

//Ke = 0.0004;
//Ke=0.001;
//Ke=0.4;%0.054;
	float Ke = 0.27;
	
	if(input > 0.9)
		current = I-I*Ke+w;
	else
		current = I -I*Ke;

	return current;
}