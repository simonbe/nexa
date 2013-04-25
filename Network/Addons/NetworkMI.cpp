#include "NetworkMI.h"
#include <math.h>


ProjectionModifierMIRateUnit::ProjectionModifierMIRateUnit(ProjectionModifier* eventHc)
{
	taupdt = 0.1;
	sprn = 1; prn = 1;
	m_eventId = 2;
	//miCi = miCj = miCij = 0;
	eventHc->AddChildProjectionModifier(this);
	m_name = "ProjectionModifierMIRateUnit";
}

ProjectionModifierMIRateUnit::~ProjectionModifierMIRateUnit()
{
	miCi.clear(); miCj.clear();
	miCij.clear();
}

void ProjectionModifierMIRateUnit::Initialize(Projection* Projection)
{
	m_projectionFixedMCs = Projection;

	for(int i = 0;i<m_parentPopulationModifier.size();i++)
	{
		m_parentPopulationModifier[i]->AddChildProjectionModifier(this);
	}

	for(int i=0;i<m_parentProjectionModifier.size();i++)
	{
		m_parentProjectionModifier[i]->AddChildProjectionModifier(this);
	}

	vector<long> localIds = m_projectionFixedMCs->GetPostLocalIds(); // minicolumn post ids
	unsigned long postNrUnits = localIds.size();//m_projectionFixed->PostLayer()->NrUnits(); // default: they are the same layer
	unsigned long preNrUnits = m_projectionFixedMCs->PreLayer()->GetNrUnitsTotal();

	miCj = vector<float>(postNrUnits,0.01);//1.0);
	miCi = vector<float>(preNrUnits,0.01);//1.0);
	miCij = vector<vector<float> >(preNrUnits,vector<float>(postNrUnits,0.01));//1.0));

	/*vector<Hypercolumn*> hsPre = ((RateUnit*)this->GetProjection()->PreLayer())->GetHypercolumns();
	vector<Hypercolumn*> hsPost = ((RateUnit*)this->GetProjection()->PostLayer())->GetHypercolumns();

	Projection* hcConn = hsPre[0]->GetProjectionFromUnitId(hsPost[0]->GetUnitId());

	if(hcConn!=NULL)
	{
		vector<ProjectionModifier*> events = hcConn->GetEvents();

		if(events.size()>0)
		{
			ProjectionModifierMIHypercolumn* e = (ProjectionModifierMIHypercolumn*)events[0]; // assuming only one type
			e->AddChildProjectionModifier(this);//AddChildEventConnMC(this);
		}
	}
	*/
}

void ProjectionModifierMIRateUnit::SetProjection(Projection* c)
{
	m_projectionFixedMCs = (ProjectionFixed*)c;
}

void ProjectionModifierMIHypercolumn::SetProjection(Projection* c)
{
	m_projectionFixedHCs = (ProjectionFixed*)c;
}

void ProjectionModifierMIRateUnit::Modify()
{
	TimingStart(m_name);

	// make all values local
	// MPILocal...
	if(IsOn() == false) return;

	if(fabs(prn)<EPS || fabs(sprn)<EPS) return;

	float prntaupdt = 0.00625/2;//taupdt*fabs(sprn);
	//prntaupdt2 = taupdt2*fabs(sprn);

	//float srctrc_i = m_pre->GetValue();
	//float srctrc_j = m_post->GetValue();
	
	vector<float> postValues = m_projectionFixedMCs->GetPostValues(); // will fetch values from units in localIds (derived from m_hasWeights, this can be optimized by a static fetch/directly from preValues, but currently most general)
	vector<float> preValues = m_projectionFixedMCs->GetPreValues((m_projectionFixedMCs->GetPostIds())[0]); // will be fully connected, so this can be independent of j
	// also, postValues and preValues will typically be the same

	for(int i=0;i<preValues.size();i++)
		miCi[i] += (preValues[i] - miCi[i])*prntaupdt; //(srctrc_i - miCi)*prntaupdt;

	for (int j=0;j<postValues.size();j++)//m_projectionFixed->PostIds()->size();j++)
	{
		miCj[j] += (postValues[j]- miCj[j])*prntaupdt;//(srctrc_j - miCj)*prntaupdt;

		for(int i=0;i<preValues.size();i++) 
		{
			miCij[i][j] += (preValues[i]*postValues[j] - miCij[i][j])*prntaupdt;
		}
	}

	TimingStop(m_name);
}

void ProjectionModifierMIHypercolumn::Modify()
{
	TimingStart(m_name);

	if(IsOn() == false) return;

	float C = 1;
	float c=C,ci,cj,cij,icij;

	//hmi = 0;
	//hji = 0;

	//for(int i=0;i<m_childProjectionModifier.size();i++)
	//{
	vector<long> postIds = m_projectionFixedHCs->GetPostIds(); // hypercolumn post ids
	
	// assuming each post has the same pres here (otherwise put in loop)
	//GetPreIdsAll // to get union
	vector<long>* preIds = m_projectionFixedHCs->GetPreIdsAll();//m_projectionFixedHCs->GetPreIds(localIds[0]);//(*m_projectionFixedHCs->PreIds())[0];
	vector<int> structurePre = m_projectionFixedHCs->PreLayer()->GetStructure();
	vector<int> structurePost = m_projectionFixedHCs->PostLayer()->GetStructure();

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

		//vector<int> mcIndexesJ = ((PopulationColumns*)m_projectionFixedHCs->PostLayer())->GetRateUnitsIndexes(localIds[j]);
		//if(j==0)		
//			startJ = mcIndexesJ[0];
		
		int currentII = 0;
		for(int i=0;i<preIds->size();i++)
		{
			//vector<int> mcIndexesI = ((PopulationColumns*)m_projectionFixedHCs->PreLayer())->GetRateUnitsIndexes(i);
			//Hypercolumn* hPre = (Hypercolumn*)(m_network->GetUnitFromId((*preIds)[i]));
			
			for(int jj=0;jj<hPost->GetRateUnits().size();jj++)//mcIndexesJ.size();jj++)
			{
				indexJ = jj+currentJJ;//mcIndexesJ[jj]-startJ;
				cj = ((ProjectionModifierMIRateUnit*)m_childProjectionModifier[0])->GetCj(indexJ);//j

				if(cj>0.0)
				{
					// currently assuming same number of minicolumns in each hypercolumn - hPre not local - but could switch to use GetStructure instead
					for(int ii=0;ii<structurePre[(*preIds)[i] - (*preIds)[0]];ii++)//hPost->GetRateUnits().size();ii++)//hPre->GetRateUnits().size();ii++)//mcIndexesI.size();ii++)
					{
						indexI = ii+currentII;//mcIndexesI[ii];
						ci = ((ProjectionModifierMIRateUnit*)m_childProjectionModifier[0])->GetCi(indexI);//i
						cij = ((ProjectionModifierMIRateUnit*)m_childProjectionModifier[0])->GetCij(indexI,indexJ);//i,j

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


vector<vector<float> > ProjectionModifierMIHypercolumn::GetValuesToRecord()
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