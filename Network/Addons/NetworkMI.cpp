#include "NetworkMI.h"
#include <math.h>


ConnectionModifierMIRateUnit::ConnectionModifierMIRateUnit(ConnectionModifier* eventHc)
{
	taupdt = 0.1;
	sprn = 1; prn = 1;
	m_eventId = 2;
	//miCi = miCj = miCij = 0;
	eventHc->AddChildConnectionModifier(this);
	m_name = "ConnectionModifierMIRateUnit";
}

ConnectionModifierMIRateUnit::~ConnectionModifierMIRateUnit()
{
	miCi.clear(); miCj.clear();
	miCij.clear();
}

void ConnectionModifierMIRateUnit::Initialize(Connection* connection)
{
	m_connectionFixedMCs = connection;

	for(int i = 0;i<m_parentPopulationModifier.size();i++)
	{
		m_parentPopulationModifier[i]->AddChildConnectionModifier(this);
	}

	for(int i=0;i<m_parentConnectionModifier.size();i++)
	{
		m_parentConnectionModifier[i]->AddChildConnectionModifier(this);
	}

	vector<long> localIds = m_connectionFixedMCs->GetPostLocalIds(); // minicolumn post ids
	unsigned long postNrUnits = localIds.size();//m_connectionFixed->PostLayer()->NrUnits(); // default: they are the same layer
	unsigned long preNrUnits = m_connectionFixedMCs->PreLayer()->GetNrUnitsTotal();

	miCj = vector<float>(postNrUnits,0.01);//1.0);
	miCi = vector<float>(preNrUnits,0.01);//1.0);
	miCij = vector<vector<float> >(preNrUnits,vector<float>(postNrUnits,0.01));//1.0));

	/*vector<Hypercolumn*> hsPre = ((RateUnit*)this->GetConnection()->PreLayer())->GetHypercolumns();
	vector<Hypercolumn*> hsPost = ((RateUnit*)this->GetConnection()->PostLayer())->GetHypercolumns();

	Connection* hcConn = hsPre[0]->GetConnectionFromUnitId(hsPost[0]->GetUnitId());

	if(hcConn!=NULL)
	{
		vector<ConnectionModifier*> events = hcConn->GetEvents();

		if(events.size()>0)
		{
			ConnectionModifierMIHypercolumn* e = (ConnectionModifierMIHypercolumn*)events[0]; // assuming only one type
			e->AddChildConnectionModifier(this);//AddChildEventConnMC(this);
		}
	}
	*/
}

void ConnectionModifierMIRateUnit::SetConnection(Connection* c)
{
	m_connectionFixedMCs = (ConnectionFixed*)c;
}

void ConnectionModifierMIHypercolumn::SetConnection(Connection* c)
{
	m_connectionFixedHCs = (ConnectionFixed*)c;
}

void ConnectionModifierMIRateUnit::Modify()
{
	TimingStart(m_name);

	// make all values local
	// MPILocal...
	if(IsOn() == false) return;

	if(prn==0 || sprn==0) return;

	float prntaupdt = 0.00625/2;//taupdt*fabs(sprn);
	//prntaupdt2 = taupdt2*fabs(sprn);

	//float srctrc_i = m_pre->GetValue();
	//float srctrc_j = m_post->GetValue();
	
	vector<float> postValues = m_connectionFixedMCs->GetPostValues(); // will fetch values from units in localIds (derived from m_hasWeights, this can be optimized by a static fetch/directly from preValues, but currently most general)
	vector<float> preValues = m_connectionFixedMCs->GetPreValues((m_connectionFixedMCs->GetPostIds())[0]); // will be fully connected, so this can be independent of j
	// also, postValues and preValues will typically be the same

	for(int i=0;i<preValues.size();i++)
		miCi[i] += (preValues[i] - miCi[i])*prntaupdt; //(srctrc_i - miCi)*prntaupdt;

	for (int j=0;j<postValues.size();j++)//m_connectionFixed->PostIds()->size();j++)
	{
		miCj[j] += (postValues[j]- miCj[j])*prntaupdt;//(srctrc_j - miCj)*prntaupdt;

		for(int i=0;i<preValues.size();i++) 
		{
			miCij[i][j] += (preValues[i]*postValues[j] - miCij[i][j])*prntaupdt;
		}
	}

	TimingStop(m_name);
}

void ConnectionModifierMIHypercolumn::Modify()
{
	TimingStart(m_name);

	if(IsOn() == false) return;

	float C = 1;
	float c=C,ci,cj,cij,icij;

	//hmi = 0;
	//hji = 0;

	//for(int i=0;i<m_childConnectionModifier.size();i++)
	//{
	vector<long> postIds = m_connectionFixedHCs->GetPostIds(); // hypercolumn post ids
	
	// assuming each post has the same pres here (otherwise put in loop)
	//GetPreIdsAll // to get union
	vector<long>* preIds = m_connectionFixedHCs->GetPreIdsAll();//m_connectionFixedHCs->GetPreIds(localIds[0]);//(*m_connectionFixedHCs->PreIds())[0];
	vector<int> structurePre = m_connectionFixedHCs->PreLayer()->GetStructure();
	vector<int> structurePost = m_connectionFixedHCs->PostLayer()->GetStructure();

	if(m_firstRun == true)
	{
		m_firstRun = false;
		hmi = vector<vector<float> >(preIds->size(),vector<float>(postIds.size(),0.0));
		hji = vector<vector<float> >(preIds->size(),vector<float>(postIds.size(),0.0));
		miDij = vector<vector<float> >(preIds->size(),vector<float>(postIds.size(),0.0));
	}

	int indexJ, indexI;
	int startJ;

	int currentJJ = 0;
	for(int j=0;j<postIds.size();j++)
	{
		Hypercolumn* hPost = (Hypercolumn*)(m_network->GetUnitFromId(postIds[j]));
		//vector<RateUnit*> mcsPost = hPost->GetRateUnits();

		//vector<int> mcIndexesJ = ((PopulationColumns*)m_connectionFixedHCs->PostLayer())->GetRateUnitsIndexes(localIds[j]);
		//if(j==0)		
//			startJ = mcIndexesJ[0];
		
		int currentII = 0;
		for(int i=0;i<preIds->size();i++)
		{
			//vector<int> mcIndexesI = ((PopulationColumns*)m_connectionFixedHCs->PreLayer())->GetRateUnitsIndexes(i);
			//Hypercolumn* hPre = (Hypercolumn*)(m_network->GetUnitFromId((*preIds)[i]));
			
			for(int jj=0;jj<hPost->GetRateUnits().size();jj++)//mcIndexesJ.size();jj++)
			{
				indexJ = jj+currentJJ;//mcIndexesJ[jj]-startJ;
				cj = ((ConnectionModifierMIRateUnit*)m_childConnectionModifier[0])->GetCj(indexJ);//j

				if(cj>0.0)
				{
					// currently assuming same number of minicolumns in each hypercolumn - hPre not local - but could switch to use GetStructure instead
					for(int ii=0;ii<structurePre[(*preIds)[i] - (*preIds)[0]];ii++)//hPost->GetRateUnits().size();ii++)//hPre->GetRateUnits().size();ii++)//mcIndexesI.size();ii++)
					{
						indexI = ii+currentII;//mcIndexesI[ii];
						ci = ((ConnectionModifierMIRateUnit*)m_childConnectionModifier[0])->GetCi(indexI);//i
						cij = ((ConnectionModifierMIRateUnit*)m_childConnectionModifier[0])->GetCij(indexI,indexJ);//i,j

						if(cij>0.0 && ci>0)
						{
							icij = log(cij*c/(ci*cj));

							// could use _finite in msc++, isinf, isnan in c99
							if(icij!=std::numeric_limits<float>::infinity() && icij != -std::numeric_limits<float>::infinity())
							{

								if(icij<MIN_FLOAT)
									icij = EPS_FLOAT;
								float Dhmi = cij/c*icij;
								float Dhji = cij/c*log(cij/c);

								hmi[i][j] += Dhmi;
								hji[i][j] -= Dhji;
							}
						}
					}
				}
			}

			if(hmi[i][j]!=std::numeric_limits<float>::infinity() && hmi[i][j] != -std::numeric_limits<float>::infinity() &&
				hji[i][j]!=std::numeric_limits<float>::infinity() && hji[i][j] != -std::numeric_limits<float>::infinity())
			{
				if(hji[i][j]>0 && hmi[i][j] > 0)
				{
					miDij[i][j] = 1 - hmi[i][j]/hji[i][j];
					//miDji = miDij // using symmetry currently not implemented
					// could also implement it by half connectivity

					if(miDij[i][j]<0.0)
						miDij[i][j] = 0.0;
				}
			}

			// assuming same number of hypercolumns - hPre not local but could switch to use GetStructure instead
			currentII+=structurePre[(*preIds)[i] - (*preIds)[0]];//hPost->GetRateUnits().size();//hPre->GetRateUnits().size();
		}

		currentJJ+=structurePost[postIds[j]];//hPost->GetRateUnits().size();
	}

	TimingStop(m_name);
}


vector<vector<float> > ConnectionModifierMIHypercolumn::GetValuesToRecord()
{
	if(IsRecording())
	{
		vector<vector<float> > data;

		for(int i=0;i<miDij.size();i++)
			data.push_back(miDij[i]);

		return data; // miDij
	}
	else return vector<vector<float> >(0);
}