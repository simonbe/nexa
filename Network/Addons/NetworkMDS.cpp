#include "NetworkMDS.h"
#include "NetworkMI.h"
#include "NetworkCorr.h"

void ProjectionModifierMDS::SetProjection(Projection* c)
{
	m_projectionFixed = (ProjectionFixed*)c;
}

void ProjectionModifierMDS::Modify()
{
	if(IsOn() == false) return;
	if(prn == 0 || sprn == 0) return;

	TimingStart(m_name);

	vector<long> localIdsPost = m_projectionFixed->GetPostLocalIds();

	LayerMDS* layerMDS = (LayerMDS*)(m_parentPopulationModifier[0]);
	m_Xi = layerMDS->GetCurrentXi();

	float dij, ddij;

	vector<vector<float> > miDij;

	if(m_usePearson == true)
	{
		ProjectionModifierPearson* e = (ProjectionModifierPearson*)m_projectionFixed->GetEvent("pearson");
		miDij = e->GetRij();
	}
	else
	{
		ProjectionModifierMIHypercolumn* e = (ProjectionModifierMIHypercolumn*)m_projectionFixed->GetEvent("mihypercolumn");
		miDij = e->GetDij();//((ProjectionModifierMIHypercolumn*)m_childProjectionModifier[0])->GetDij();//m_MI->GetDij();
	}

	float meanddij = 0;
	Xm = vector<float>(m_mdsDim);

	for(int d=0;d<m_mdsDim; d++)
	{
		//oldXi_i[d] = Xi_i[d];
		Xm[d] = 0.0; // for re-position
	}

	vector<vector<float> > diffXi = vector<vector<float> >(miDij.size(),vector<float>(m_mdsDim));

//	bool useThreshold = false;
	// assuming here that calculation of the matrix NxN is divided up like mxN
	float totLocalStress = 0;

	for(int i=0;i<miDij[0].size();i++)
	{
		int iIndex = (int)localIdsPost[i];

		for(int j=0;j<miDij.size();j++)
		{
			dij = 0;

			// euclidean distance
			for(int d=0;d<m_mdsDim; d++)
				dij+=((*m_Xi)[j][d] - (*m_Xi)[iIndex][d])*((*m_Xi)[j][d] - (*m_Xi)[iIndex][d]);//dij+=pow((*m_Xi)[j][d] - (*m_Xi)[iIndex][d],2);//dij += pow(Xi_j[d] - Xi_i[d],2);

			dij = sqrt(dij);
			totLocalStress+=dij;

			bool doChange = true;

			if(m_useThreshold)
			{
				//float threshold = 0.9;
				
				if(miDij[j][i] > m_threshold)
				{
					if(dij>m_threshold)
					{
						doChange = false;
					}
				}
			}

			if(doChange == true)
			{
				// difference to real value
				//ddij = (dij-miDij[j][i]);
				ddij = (dij-miDij[j][i]);///(miDij[j][i]+0.1); // avoid division by zero
				//ddij = fabs((dij-miDij[j][i])/(miDij[j][i]+1));
				meanddij += ddij*ddij;//pow(ddij,2);
				ddij *= m_MDSK/10;
				//ddij /= 2/(m_MDSK/10);
				//ddij = 0.001;

				for(int d=0;d<m_mdsDim;d++)
				{
					//diffXi[iIndex][d] += ddij*((*m_Xi)[j][d]-(*m_Xi)[iIndex][d]);	//Xi_i[d] += ddij*(Xi_j[d]-Xi_i[d]);
					//diffXi[j][d] -= ddij*((*m_Xi)[j][d]-(*m_Xi)[iIndex][d]);
					float v=ddij*((*m_Xi)[j][d]-(*m_Xi)[iIndex][d]);
					diffXi[iIndex][d] += v;	//Xi_i[d] += ddij*(Xi_j[d]-Xi_i[d]);
					diffXi[j][d] -= v;

					// if...
					//	Xm[d]+= Xi_i[d]; // only be done on one Projection from hc
				}
			}
		}
	}

	layerMDS->AddDiffXi(diffXi);
	layerMDS->SetLocalStress(totLocalStress);

	TimingStop(m_name);
}

void ProjectionModifierMDS::TranslateToOrigin(vector<float> XmTot)
{
	for(int d=0;d<m_mdsDim;d++) 
	{
		Xi_i[d] -= XmTot[d];
		Xi_j[d] -= XmTot[d];

		meands += (Xi_i[d]-oldXi_i[d])*(Xi_i[d]-oldXi_i[d]);//pow(Xi_i[d]-oldXi_i[d],2);
	}

	meands = 0;
}

void ProjectionModifierMDS::AddParentPopulationModifier(PopulationModifier* e)
{
	m_parentPopulationModifier.push_back(e);
	m_mdsDim = ((LayerMDS*)e)->GetDimension();
	Xi_i = vector<float>(m_mdsDim,0.0);
	Xi_j = vector<float>(m_mdsDim,0.0);
	for(int i=0;i<Xi_i.size();i++)
	{
		Xi_i[i] = (rand()/(float)RAND_MAX) - 0.5;
		Xi_j[i] = (rand()/(float)RAND_MAX) - 0.5;
	}

	oldXi_i = vector<float>(m_mdsDim,0.0);
}

void LayerMDS::Simulate()
{
	TimingStart(m_name);

	if(IsOn() == false) return;

	// sum all tot diffs
	TimingStart("CommunicationLayerMDS");
	m_diffXi = ((PopulationColumns*)m_population)->MPI()->MPILayerReduceVariable(m_diffXi);
	m_totStress = (((PopulationColumns*)m_population)->MPI()->MPILayerReduceVariable(vector<float>(1,m_localStress)))[0];
	TimingStop("CommunicationLayerMDS");

	// replace Xi
	for(int i=0;i<m_Xi.size();i++)
	{
		for(int j=0;j<m_Xi[0].size();j++)
		{
			m_Xi[i][j] += m_diffXi[i][j];
		}
	}

	
	m_totChange = 0;
	for(int i=0;i<m_Xi.size();i++)
	{
		for(int j=0;j<m_Xi[0].size();j++)
		{
			m_totChange+=fabs(m_diffXi[i][j]);
		}
	}
	if(m_population->network()->MPIGetNodeId() == 0 && DEBUG_LEVEL > 2)
	{
		cout<<"MDS-totDiff="<<m_totChange<<" ";
	}
	if(fabs(m_totChange-m_prevDist)<MDS_EPS) {
		m_distUnchanged++;
	}
	else {
		m_distUnchanged=0;
		m_prevDist=m_totChange;
	}
	SwitchOnOff(m_distUnchanged<MDS_NUM);

	// reset
	m_diffXi = vector<vector<float> >(m_diffXi.size(),vector<float>(m_diffXi[0].size()));

	if(this->network()->MPIGetNodeId() == 0)
	{
		if(IsRecording("stress")) // how slow? can put in hash to speed up
		{
			vector<float> v(2);
			v[0] = this->network()->GetCurrentTimeStep();
			v[1] = m_totChange;

			m_recordedValues.push_back(v);
		}
	}

	TimingStop(m_name);
}

void ProjectionModifierMDS::Dispose()
{
}