#include "NetworkSanger.h"

ProjectionModifierSanger::ProjectionModifierSanger()
{
	m_eventId = 14;
	m_etaHebb = 0.005;
	m_transferFunction = new TransferLinear(false);
}

ProjectionModifierSanger::ProjectionModifierSanger(float etaHebb)
{
	m_eventId = 14;
	m_etaHebb = etaHebb;
	m_transferFunction = new TransferLinear(false);
}

/*void ProjectionModifierSanger::SetEtaHebb(float etaHebb)
{
	m_etaHebb = etaHebb;
}*/

void ProjectionModifierSanger::Initialize(Projection* Projection)
{
	network(Projection->network());

	m_projectionFixed = Projection;
	m_firstRun = true;
	m_idsPost = m_projectionFixed->GetPostIds();
}

void ProjectionModifierSanger::SetProjection(Projection* c)
{
	m_projectionFixed = (ProjectionFixed*)c;
}

void ProjectionModifierSanger::Modify()
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

	for(int j=0;j<postValues.size();j++)
	{
		vector<float> preValues = m_projectionFixed->GetPreValues(m_idsPost[j]);
		postId = m_idsPost[j];
		y = postValues[j];

		for(int i=0;i<preValues.size();i++)
		{
			preId = (*preIds)[j][i];
			x = preValues[i];

			if(preId == postId)
				bool b = false;

			weight = network()->GetWeight(preId,postId);
			weight = weight + m_etaHebb*y*(x-y*weight); // oja's, only first pc
			//weight - phi*y*weight;//w[count]=w[count]-phi[n]*y[n]*w[count]

			network()->SetWeight(weight,preId,postId);
		}

		// weight normalization

		/*totWeights = 0;

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
		}*/
	}
}

void ProjectionModifierSanger::Simulate(UnitModifier* e)
{

}

/*void ProjectionModifierTriesch::Clear()
{
	// 1. Clear weight values

	vector<float> postValues = m_projectionFixed->GetPostValues();
	vector<vector<long> >* preIds = m_projectionFixed->PreIds(); // move to initializer (until Invalidate().. called)

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
}*/

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
        elif learning_rule==3: # hebb
            count=0
            for n from 0<=n < num_neurons:
                phi[n]=eta*y[n]
                
                for i from 0<=i<num_inputs:
                    w[count]=w[count]+phi[n]*x[count]-eta*decay*w[count]
                    count=count+1
                    
                th[n]=th[n]+(y[n]*y[n]-th[n])/tau
        elif learning_rule==4: # user1  
            #  start from equation 14, replace Y with a+bV+c*V**2, where V is the
            #postsynaptic activity
            #   S is sigmoid of input
            

            count=0
            for n from 0<=n < num_neurons:
                for i from 0<=i<num_inputs:
                    tmp1=a+b*y[n]+c*y[n]*y[n]
                    
                    w[count]=w[count]+eta*( (theta_o/th[n])**3 * 
                                        (x[count]*tmp1)**2 *
                                        (x[count]*tmp1-th[n])
                                      
                                        -(theta_o/theta[n])*
                                          (x[count]*tmp1)*w[count])-eta*decay*w[count]
                    
                    count=count+1
                    
                th[n]=th[n]+(y[n]*y[n]-th[n])/tau
        elif learning_rule==5: # user2: hebb with scaling
            count=0
            for n from 0<=n < num_neurons:
                phi[n]=eta*y[n]
                
                for i from 0<=i<num_inputs:
                    w[count]=w[count]+(phi[n]*x[count]-eta*decay*w[count])/th[n]
                    count=count+1
                    
                th[n]=th[n]+(y[n]*y[n]-th[n])/tau
        elif learning_rule==6: # unified 1
            count=0
            for n from 0<=n < num_neurons:
                for i from 0<=i<num_inputs:
                    w[count]=w[count]+ (
                       eta*(theta_o**theta_p/th[n]**theta_p)*x[count]*
                               (x[count]*y[n])**2*(x[count]*y[n]-th[n])-
                       unified_decay*eta*theta_o/th[n]*x[count]**2*y[n]*w[count]
                    )
                    
                    count=count+1
                    
                th[n]=th[n]+(y[n]*y[n]-th[n])/tau
            
            
            
        elif learning_rule==7: # unified 2
            count=0
            for n from 0<=n < num_neurons:
                for i from 0<=i<num_inputs:
                    w[count]=w[count]+ (
                       eta*(theta_o**theta_p/th[n]**theta_p)*x[count]*
                               (x[count]*y[n])**2*(y[n]-th[n])-
                       unified_decay*eta*theta_o/th[n]*x[count]*y[n]*w[count]
                    )
                    
                    count=count+1
                    
                th[n]=th[n]+(y[n]*y[n]-th[n])/tau
            
            
            
        else:
            raise ValueError
        
                
        if weight_stabilization_type==0: # nothing
            pass
        elif weight_stabilization_type==1: # oja
            count=0
            for n from 0<=n < num_neurons:
                for i from 0<=i<num_inputs:
                    w[count]=w[count]-phi[n]*y[n]*w[count]
                    count=count+1
        elif weight_stabilization_type==2: # strict
            count=0
            for n from 0<=n < num_neurons:
                sum=0.0
                for i from 0<=i<num_inputs:
                    sum=sum+w[count]*w[count]
                    count=count+1
                    
                count=count-num_inputs
                sum=sum/sqrt(sum)
                
                for i from 0<=i<num_inputs:
                    w[count]=w[count]/sum
                    count=count+1
            
        elif (weight_stabilization_type==3 or weight_stabilization_type==6): # saturation
            count=0
            for n from 0<=n < num_neurons:
                for i from 0<=i<num_inputs:
                    if w[count]>weight_stabilization_top:
                        w[count]=weight_stabilization_top
                    elif w[count]<weight_stabilization_bottom:
                        w[count]=weight_stabilization_bottom
                    count=count+1
        elif weight_stabilization_type==4: # weight decay
            pass  # done above
        elif weight_stabilization_type==5: # saturation w/o zero cross
            pass  # not implemented
        elif weight_stabilization_type==6: # weight decay with saturation
            pass # done above
        else: 
            raise ValueError*/