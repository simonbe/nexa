#include "NetworkKussul.h"

ConnectionModifierKussul::ConnectionModifierKussul()
{
	m_eventId = 7;
	m_transferFunction = new TransferLinear(false);
}

void ConnectionModifierKussul::Simulate(UnitModifier* e)
{

}

void ConnectionModifierKussul::Modify()
{
	if(IsOn() == false) return;

	if(m_connectionFixed->PreLayer()->network()->MPIGetNodeId() == 0) // will be put in mpi-class
	{
		cout<<".";
		cout.flush();
	}

	PopulationColumns* layer = (PopulationColumns*)m_connectionFixed->PostLayer();//m_population;
	
	layer->MPI()->MPIMakeHypercolumnsValuesLocal(); // distribute all post layer values and keep
	vector<RateUnit*> minicolumns = ((PopulationColumns*)m_connectionFixed->PostLayer())->GetRateUnits();
	vector<float> supOutput(minicolumns.size());
	for(int i=0;i<minicolumns.size();i++)
		supOutput[i] = minicolumns[i]->GetValue();

	vector<float> postValues = m_connectionFixed->GetPostValues();
	
	// run forward
	vector<float> output(postValues.size());
	vector<vector<long> >* preIds;// = m_connectionFixed->PreIds();

	m_idsPost = m_connectionFixed->GetPostIds();

	for(int j=0;j<postValues.size();j++)
	{
		vector<float> preValues = m_connectionFixed->GetPreValues(m_idsPost[j]);

		for(int i=0;i<preValues.size();i++)
		{	
			long preId = (*preIds)[j][i];
			float weight = network()->GetWeight(preId,m_idsPost[j]);

			output[j] += preValues[i] * weight;
		}
	}

	for(int i=0;i<output.size();i++)
	{
		((RateUnit*)m_connectionFixed->PostLayer()->network()->GetUnitFromId(m_idsPost[i]))->SetValue(output[i]);
	}

	layer->MPI()->MPIMakeHypercolumnsValuesLocal(); // MPI distributor

	int wIndex = -1;
	float wValue = -10e8;
	int swIndex = -1;
	float swValue = -10e8;

	// determine winner, v
	minicolumns = ((PopulationColumns*)m_connectionFixed->PostLayer())->GetRateUnits();
	
	for(int i=0;i<minicolumns.size();i++)
	{
		if(minicolumns[i]->GetValue()>wValue)
		{
			wValue = minicolumns[i]->GetValue();
			wIndex = i;
		}

		if(supOutput[i]>swValue)
		{
			swValue = supOutput[i];
			swIndex = i;
		}
	}

	// if v==c, do nothing
	if(wIndex == swIndex)
	{
		bool b=false;
		//nrCorrect++;
	}
	else // if v!=c
	{
		// w_ic = w_ic + a_i
		if(minicolumns[swIndex]->IsLocal() == true)
		{
//			int j = minicolumns[swIndex]->GetUnitIdLocal();
			long unitId = minicolumns[swIndex]->GetUnitId();

			// will be deep copies, change?
			vector<long> preIds = m_connectionFixed->GetPreIds(unitId);//PreIds();
			vector<float> preValues = m_connectionFixed->GetPreValues(unitId);

			for(int i=0;i<preIds.size();i++)
			{	
				if(preValues[i]>0.5)
				{
					long preId = preIds[i];//(*preIds)[j][i];
					float weight = network()->GetWeight(preId,unitId);//m_idsPost[j]);
					weight++;
					if(weight>1)
						bool dbg = true;
					network()->SetWeight(weight,preId,unitId);//m_idsPost[j]);
				}
			}
		}
		
		// w_iv = w_iv - a_i
		if(minicolumns[wIndex]->IsLocal() == true)
		{
			long unitId = minicolumns[wIndex]->GetUnitId();

			vector<long> preIds = m_connectionFixed->GetPreIds(unitId);//PreIds();
			vector<float> preValues = m_connectionFixed->GetPreValues(unitId);

			for(int i=0;i<preIds.size();i++)
			{	
				if(preValues[i]>0.5)
				{
					long preId = preIds[i];//(*preIds)[j][i];
					float weight = network()->GetWeight(preId,unitId);//m_idsPost[j]);
					weight--;
					if(weight>0)
						network()->SetWeight(weight,preId,unitId);//m_idsPost[j]);
					else
						network()->SetWeight(0.0,preId,unitId);//m_idsPost[j]);
				}
			}
		}

		// where a_i output signal (0 or 1)

		// if w_iv < 0, w_iv = 0;
	}
}


void ConnectionModifierKussul::Initialize(Connection* connection)
{
	network(connection->network());

	m_connectionFixed = connection;
	vector<float> preValues = vector<float>(m_connectionFixed->PreLayer()->GetUnits().size(),0.0);//m_connectionFixed->GetPreValues();
	vector<float> postValues = m_connectionFixed->GetPostValues();

	//m_Ai = vector<float>(preValues.size(),0.01);//= m_Aj = 0.001;
	//m_Aj = vector<float>(postValues.size(),0.01);
	//m_Aij = vector<vector<float> >(preValues.size(),vector<float>(postValues.size(),0.01*0.01));//0.001*0.001;
	//m_beta = vector<float>(postValues.size(),0);//0;

	m_firstRun = true;
	m_idsPost = m_connectionFixed->GetPostIds();
}

void ConnectionModifierKussul::SetConnection(Connection* c)
{
	m_connectionFixed = (ConnectionFixed*)c;
}