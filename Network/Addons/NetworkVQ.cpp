#include <ctime>
#include "NetworkVQ.h"
#include "NetworkMDS.h"

/// puts the code vectors at random data points
void LayerVQ::InitiateVQUnits()
{
	vector<int> selected;
	vector<int> indexes(m_data->size());
	m_data = m_childPopulationModifier[0]->GetOutput(); // latest copy of (mds) data

	for(int i=0;i<m_data->size();i++)
	{
		indexes[i] = i;
	}

	unsigned int seed;
	if(m_population->network()->MPIGetNodeId() == 0)
	{
		seed = (unsigned)time(NULL); // only set globally, may mess up
		MPI_Bcast(&seed,1,MPI_UNSIGNED,0,NETWORK_COMM_WORLD);
	}
	else
	{
		MPI_Bcast(&seed,1,MPI_UNSIGNED,0,NETWORK_COMM_WORLD);
	}

	//srand(seed);

	for(int i=0;i<trgnhyp;i++)
	{
		int r = rand() % indexes.size();
		selected.insert(selected.end(),indexes[r]);
		indexes.erase(indexes.begin() + r);	
	}

 	// select as initial code vectors
	for(int i=0;i<selected.size();i++)
	{
		for(int j=0;j<m_MDSDIM;j++)
		{
			vqu[i][j] = (*m_data)[selected[i]][j];
		}
	}
}

void LayerVQ::Simulate()
{
	TimingStart(m_name);

	if(IsOn() == false) return;

	// point to latest copy
	m_data = m_childPopulationModifier[0]->GetOutput();

	if(m_firstRun == true)
	{
		m_mpiSize = m_population->network()->MPIGetNrProcs();
		m_mpiRank = m_population->network()->MPIGetNodeId();
	}

	if(m_type == LayerVQ::VQStandard)
	{
		if(m_firstRun == true)
		{
			m_mpiNrNodes = m_population->network()->MPIGetNrProcs();
			srcnhyp = m_data->size();

			// retrieve parts
			vector<vector<long> > parts = m_population->MPI()->DivideEqualByUnits(srcnhyp,m_population->network()->MPIGetNodeId(),m_population->network()->MPIGetNrProcs());
			m_mpiParts = parts[m_population->network()->MPIGetNodeId()];
			m_sizeThisNode = m_mpiParts[1] - m_mpiParts[0];

			m_MDSDIM = (*m_data)[0].size();
			vqwin = vector<int>(srcnhyp);
			vqu = vector<vector<float> >(trgnhyp,vector<float>(m_MDSDIM));
			InitiateVQUnits();
			vquDiff = vqu;
			vqd = vector<float>(trgnhyp);
			vqn = vqd;
			m_distFracChanged = false;
			m_distFrac = 0;
		}

		simstep = m_population->network()->GetCurrentTimeStep();

		updvq();
	}
	else if(m_type == LayerVQ::VQCSL)
	{
		if(m_firstRun == true)
		{
			//m_csl = new CSL();
			m_csl->SetMPIParameters(this->network()->MPIGetNodeId(),this->network()->MPIGetNrProcs());
			m_csl->Initialize(m_data,m_nrCodeVectors, this->network()->MPIGetNodeId(),this->network()->MPIGetNrProcs(), true, this, NETWORK_COMM_WORLD); // all nodes involved
			//CSL_Initiate();
		}
		SwitchOnOff(m_csl->Step(m_data));
		if(this->m_nrOverlaps > 1)
			m_totalWinnersUnits = m_csl->GetTotalWinnersUnits(m_nrOverlaps);
		else
			m_totalWinnersUnits = m_csl->GetTotalWinnersUnits(); // change to pointer
		//CSL_Step();

		// order them so that values will be retrieved ordered
		for(int i=0;i<m_totalWinnersUnits.size();i++)
		{
			sort(m_totalWinnersUnits[i].begin(),m_totalWinnersUnits[i].end());
		}
		


		if(this->network()->MPIGetNodeId() == 0)
		{
			if(IsRecording("stress")) // how slow? can put in hash to speed up
			{
				//		vector<float> values(1);
				//		values[0] = m_storedValues;//S_.y_[V_M];
				vector<float> v(2);
				v[0] = this->network()->GetCurrentTimeStep();
				v[1] = m_csl->GetDistortion();

				m_recordedValues.push_back(v);
			}
		}
	}

	if(m_firstRun == true)
		m_firstRun = false;

	TimingStop(m_name);
}

float LayerVQ::euclid(vector<float> v1,vector<float> v2,int n) 
{
	float sqsum = 0;
	for (int i=0; i<n; i++) sqsum += (v1[i]-v2[i])*(v1[i]-v2[i]);
	return sqrt(sqsum);
}

void LayerVQ::updvqwin(int u)
{
	//if (prn==0 || sprn==0) return;
	// Update distances to vqu:s
	vqdst = vector<float>(trgnhyp);

	for (int v=0; v<trgnhyp; v++)
		vqdst[v] = euclid(vqu[v],(*m_data)[u],m_MDSDIM);
	// Find NWIN vq winners
	float min;
	int nw = 0;
	for (int w=0; w<m_NWIN; w++) {
		min = 1e30f;
		for (int v=0; v<trgnhyp; v++)
			if (0<=vqdst[v] && vqdst[v]<min) {
				vqwin[nw] = v;
				min = vqdst[v];
			}
		vqdst[vqwin[nw]] = -vqdst[vqwin[nw]];
		nw++;
	}
	for (int v=0; v<trgnhyp; v++) if (vqdst[v]<0) vqdst[v] = -vqdst[v];
}

float r01() {
	return (rand()/(float)RAND_MAX);
}

void LayerVQ::updvq() {

	if (fabs(prn)<EPS || fabs(sprn)<EPS) return;
	float dmin,dmax,di;
	int vmin,vmax;
	for (int v=0; v<trgnhyp; v++) { vqd[v] = 0; vqn[v] = 0; }
	/*for (int h=0; h<srcnhyp; h++)
		for (int v=0; v<trgnhyp; v++) {
			oldhuse[h][v] = huse[h][v];
			huse[h][v] = 0;
		}*/

	for(int i=0;i<vquDiff.size();i++)
		for(int j=0;j<vquDiff[0].size();j++)
			vquDiff[i][j] = 0.0;

	m_winnersUnits.clear();
	m_winnersUnits = vector<vector<int> >(trgnhyp, vector<int>(0));

	for(int h=m_mpiParts[0];h<m_mpiParts[1];h++)//for (int h=0; h<srcnhyp; h++) {
	{
		updvqwin(h);
		m_winnersUnits[vqwin[0]].push_back(h);

		/* Move the winning unit, allocate the others */
		int iwin;
		for (int w=0; w<m_NWIN; w++) {
			iwin = vqwin[w];
			if (w==0) 
			{
				for (int d=0; d<m_MDSDIM; d++)
				{
					//vqu[iwin][d] += m_MDSK * ((*m_data)[h][d] - vqu[iwin][d]);
					vquDiff[iwin][d] += m_MDSK/(float)(m_mpiNrNodes) * ((*m_data)[h][d] - vqu[iwin][d]);// nr nodes involved
				}

				vqd[iwin] += vqdst[iwin];
				vqn[iwin]++;
				//huse[h][iwin] = 2;
			} //else huse[h][iwin] = 1;
		}
	}

	// bring together
	// sum all vquDiff
	vquDiff = m_population->MPI()->MPILayerReduceVariable(vquDiff);

	// add to vqu
	for(int i=0;i<vqu.size();i++)
		for(int j=0;j<vqu[0].size();j++)
		{
			vqu[i][j]+=vquDiff[i][j];
		}

	// sum all vqd and vqn
	vqd = m_population->MPI()->MPILayerReduceVariable(vqd);
	vqn = m_population->MPI()->MPILayerReduceVariable(vqn);
	
	/* Check if some huse changed and update accordingly */
/*	for (int hi=0; hi<srcnhyp; hi++) {
		for (int hj=0; hj<trgnhyp; hj++) {
			oldhuse[hi][hj] = (huse[hi][hj] - oldhuse[hi][hj]);
			if (oldhuse[hi][hj]<0) initc(hi,hj);
		}
	}
	*/

	// Same on all nodes - could be changed

	for (int v=0; v<trgnhyp; v++) if (vqn[v]>0) vqd[v] /= vqn[v];


	dmin = vqd[0]; vmin = 0;
	dmax = vqd[0]; vmax = 0;
	for (int v=1; v<trgnhyp; v++) {
		if (vqd[v]<dmin) { dmin = vqd[v]; vmin = v; } else
			if (vqd[v]>dmax) { dmax = vqd[v]; vmax = v; }
	}

	if(m_population->network()->MPIGetNodeId() == 0 && DEBUG_LEVEL > 2) cout << "VQ-dmax="<<dmax << " ";
	int nx = (int)(3/sprn);//(int)(100/sprn);

	if ((simstep+1)%nx==0) 
	{
		/* Move unit with lowest distorsion to duplicate unit with highest resolution */

		if(!(fabs(dmin)<EPS) && m_distFracChanged == false)
			m_distFracChanged = true;

		if((m_distFrac - dmin/dmax)/m_distFrac > 1e-8 || m_distFracChanged == false) // unless no more change
		{
			//		if (dmax < pow(9,1.0f/MDSDIM)*dmin) return;
			if (sprn>0 && dmin < m_RLIM*dmax) {
				//cerr << vmin << " " << vmax << " ";
				for (int d=0; d<m_MDSDIM; d++) {
					di = 1e-6f*((float)r01() - 0.5f);
					vqu[vmin][d] = vqu[vmax][d] + di;
					vqu[vmax][d] = vqu[vmax][d] - di;
				}
			}

			m_distFrac = dmin/dmax;
		}
	}


}



///////////////////////////////////////////////////////////////////////////////
/////// CSL Implementation
///////////////////////////////////////////////////////////////////////////////

void LayerVQ::CSL_Initiate()
{
	m_x = *m_childPopulationModifier[0]->GetOutput(); // latest copy of (mds) data

	m_epsilon = 0.01;//0.001;
	m_eta = 0.01;//0.0001;//0.5;//0.001;
	
	// set initial values
	m_M = 1;
	m_currentDistortion = 1e10;
	int mpiRank = this->m_population->network()->MPIGetNodeId();
	
	// randomly select nrNeurons different code vectors

	vector<int> selected;
	vector<int> indexes(m_x.size());

	for(int i=0;i<m_x.size();i++)
	{
		indexes[i] = i;
	}

	unsigned int seed;
	if(mpiRank == 0)
	{
		seed = 1;//(unsigned)time(NULL);
		MPI_Bcast(&seed,1,MPI_UNSIGNED,0,NETWORK_COMM_WORLD); // should not be world
	}
	else
	{
		MPI_Bcast(&seed,1,MPI_UNSIGNED,0,NETWORK_COMM_WORLD);
	}

	//srand(seed);//(unsigned)time(NULL));
	// randomization also in moving the codevectors (selection step)

	if(mpiRank == 0)
	{
		cout<<mpiRank<<": m_x.size()="<<m_x.size()<<"\n";
		cout<<mpiRank<<": indexes.size()="<<indexes.size()<<"\n";
		cout<<mpiRank<<": m_nrCodeVectors="<<m_nrCodeVectors<<"\n";
		cout.flush();
	}

	for(int i=0;i<m_nrCodeVectors;i++)
	{
		// select indexes
		int r = rand() % indexes.size();
		
		selected.insert(selected.end(),indexes[r]);
		indexes.erase(indexes.begin() + r);	
	}

 	// select as initial code vectors
	for(int i=0;i<selected.size();i++)
	{
		m_c.insert(m_c.end(),m_x[selected[i]]);
	}

	if(mpiRank == 0)
	{
		cout<<"CSL_Initiate complete.\n";
	}

	// Selection sizes (declining)
	for(int i=m_c.size()-m_c.size()/2; i>0; i--)
	{
		for(int j=0;j<10;j++)
		{
			m_S.insert(m_S.end(),i);
			if(m_S.size()<m_x.size());
		}
	}

	m_M = 0; // Selection time step index
}

vector<float> LayerVQ::CSL_NewCodeVector(vector<float> oldCodeVector, vector<float> inputVector)
{
  vector<float> newVector(oldCodeVector.size());
  
  //float eta = 0.01;
  for(int i=0;i<oldCodeVector.size();i++) {
    newVector[i] = oldCodeVector[i] + m_eta*(inputVector[i] - 
					   oldCodeVector[i]);
  }

  return newVector;
}


void LayerVQ::CSL_Selection(vector<int> indexes, vector<float> D)
{
  int s = indexes.size();
  
  // 1. Compute fitness measure using subdistortion
  float gamma = 0.5;

  vector<float> g(s);

  float sum = 0;
  for(int i=0;i<s;i++)
  {
	  sum += pow(D[indexes[i]],gamma);
  }

  for(int i=0;i<s;i++)
  {
	  g[i] = pow(D[indexes[i]],gamma) / sum;
  }
  
  // 2. Determine nr of neurons to be reproduced
  // i) compute integer part

  vector<int> u(s);
  for(int i=0;i<s;i++)
  {
	  u[i] = (int)(g[i] * (float)s);
  }

  // ii) for top s - sum(u) wrt gs - u , add one
  
  int usum = 0;
  for(int i=0;i<u.size();i++) {
    usum += u[i];
  }

  vector<float> data(u.size());
  for(int i=0;i<u.size();i++) {
    data[i] = g[i]*(float)s - (float)u[i];
  }

  int nrValues = s - usum;

  vector<int> top = CSL_GetHighestIndexes(false, nrValues,data);
  
  // add one
  for(int i=0;i<top.size();i++) {
    u[top[i]] = u[top[i]] + 1;
  }

  // 3. For every j, reproduce u_j neurons by adding random perturbations (and remove if u_j == 0)
  vector<vector<float> > newC;

  /*
  unsigned int seed = GetSeed();
  
  srand(seed);//(unsigned)time(NULL));
  */
  float epsilonMove = 0.01;

  for(int i=0;i<u.size();i++) 
  {
	//  m_c.erase(indexes[i]);

	  for(int j=0;j<u[i];j++) 
	  {
		  vector<float> c(m_c[0].size());

		  for(int x=0;x<m_c[0].size();x++) 
		  {
			  float dRand = ((float)rand()/(float)RAND_MAX);

			  c[x] = (1+epsilonMove*(dRand - 0.5))*m_c[indexes[i]][x];
		  }

		  newC.insert(newC.end(),c);
	  }
  }

  //TRACK: newC empty which raises exceptions from this Loop
  for(int i=0;i<indexes.size();i++) 
  {
	  for(int j=0;j<m_c[0].size();j++)
		m_c[indexes[i]][j] = newC[i][j];
  }
}

vector<int> LayerVQ::CSL_GetHighestIndexes(bool highestAndLowest, int nrValues, vector<float> data)
{
	vector<float> sorted;
	vector<int> indexesSorted;

	// sort
	for(int i=0;i<data.size();i++) 
	{
		bool placed = false;
		for(int j=0;j<sorted.size();j++) 
		{
			if(data[i]>sorted[j]) 
			{
				sorted.insert(sorted.begin() + j,data[i]);
				indexesSorted.insert(indexesSorted.begin() + j,i);
				placed = true;
				break;
			}
		}

		if(placed==false) {
			sorted.insert(sorted.end(),data[i]);
			indexesSorted.insert(indexesSorted.end(),i);
		}
	}

	// get nrValues first ones

	//Get the highest values
	if(highestAndLowest == false) {
		vector<int> out(nrValues);

		for(int i=0;i<nrValues;i++) {
			out[i] = indexesSorted[i];
		}

		return out;
	} else {
		// get alternatively highest and lowest values
		vector<int> out(nrValues);

		int indexHigh = 0;
		int indexLow = data.size()-1;//nrValues-1;

		for(int i=0;i<nrValues;i++) {
			if(i%2 == 0) { //even 
				out[i] = indexesSorted[indexHigh];
				indexHigh++;
			} else {
				out[i] = indexesSorted[indexLow];
				indexLow--;
			}
		}

		return out;
	}
}

unsigned int LayerVQ::GetSeed()
{
	bool VQSeedUseTime = false;
	if(VQSeedUseTime == true)
	{
		int seed;

		if(this->m_population->network()->MPIGetNodeId() == 0)
		{
			seed = (unsigned)time(NULL);
			MPI_Bcast(&seed,1,MPI_UNSIGNED,0,NETWORK_COMM_WORLD);
		}
		else
		{
			MPI_Bcast(&seed,1,MPI_UNSIGNED,0,NETWORK_COMM_WORLD);
		}

		return seed;
	}
	else
	{
		return 0;
	}
}

void LayerVQ::CSL_Step()
{
	m_x = *m_childPopulationModifier[0]->GetOutput(); // latest copy of (mds) data

	m_distortion = 0;
	vector<float> D(m_c.size());
	
	//adjust vectors using competative learning

	CSL_ParallelCompetitiveLearning();
	CSL_ParallelDistortionCalculation( &m_distortion, &D );

	//3. Selection
	if(fabs(m_distortion) < EPS) // usually never happens
	{
		if(this->m_population->network()->MPIGetNodeId() == 0)
		{
			cout<<"\ndistortion == 0!";
			cout.flush();
		}

		//converged = true;
		m_currentDistortion = m_distortion;
	}
	else
	{
		if(m_S[m_M]>0)
		{
			//get the m_s[m_M] highest/lowest indexes wrt subdistortion

			vector<int> highestDIndexes = CSL_GetHighestIndexes(true,m_S[m_M],D);
			CSL_Selection(highestDIndexes, D);

			//if(count % m_settings.cslpCount == 0)
			//{
				m_M+=1;//2
				//if(m_M<0)
				//	m_M = 0;
				//if(m_M>m_S.size()-1)
				//	m_M = m_S.size()-3;
			//}
		}
	}

	m_currentDistortion = m_distortion; // not needed to be global atm.
	if(m_mpiRank == 0 && DEBUG_LEVEL > 2) 
	{
		cout << "VQ-dtot="<<m_distortion<< "(m_S[m_M]="<<m_S[m_M]<<") ";
	}

	// for connectivity
	m_totalWinnersUnits.clear();
	m_totalWinnersUnits = vector<vector<int> >(m_c.size(),vector<int>(0));

	for(int i=0;i<m_total_winners.size();i++)
	{
		m_totalWinnersUnits[m_total_winners[i]].push_back(i);
	}
}

float LayerVQ::CSL_RRValue(const vector<float>& x1, const vector<float>& x2)
{
  float val = 0;
  
  for(int i=0;i<x1.size();i++) {
    val+=(x1[i]-x2[i])*(x1[i]-x2[i]);//pow(x1[i]-x2[i],2);
  }

  return sqrt(val);
}


vector<float> LayerVQ::CSL_GetClosestCodeVectorIndexAndValue(const vector<float>& data, const vector<vector<float> >& codeVectors)
{
  vector<float> out;
  
  float minValue = 1e10;
  float minIndex = -1;
  
  for(int i=0;i<codeVectors.size();i++) {
    float val = CSL_RRValue(data,codeVectors[i]);//m_x[m],m_c[n]);
    if(val<minValue) {
      minValue = val;
      minIndex = i;
    }
  }

  out.insert(out.end(),minIndex);
  out.insert(out.end(),minValue);
  
  return out;
}


void LayerVQ::CSL_ParallelDistortionCalculation(float *distortion, vector<float> *D)
{
	int n_s = m_x.size()/m_mpiSize;

	vector<float> totD(m_c.size());
	vector<float> subD(m_c.size());

	float newDistortion = 0;
	
	int max = (m_mpiRank*n_s)+n_s;
	if(m_mpiRank == m_mpiSize-1)
		max = m_x.size();

	for(int i = m_mpiRank*n_s; i < max; i++) 
	{
		vector<float> v =  CSL_GetClosestCodeVectorIndexAndValue(m_x[i], m_c);
		subD[(int)v[0]] += v[1];
		newDistortion += v[1];
	}

 	MPI_Allreduce(&newDistortion, distortion, 1, MPI_FLOAT, MPI_SUM, NETWORK_COMM_WORLD);
	MPI_Allreduce(&subD[0], &totD[0], m_c.size(), MPI_FLOAT, MPI_SUM, NETWORK_COMM_WORLD);
	
	for(int i=0;i<totD.size();i++)
		totD[i] = totD[i]/(float)m_x.size();

	*D = totD;
	*distortion = *distortion / (float)m_x.size();
	//MPI_Barrier(NETWORK_COMM_WORLD);
	//std::cout << "Rank " << m_mpiRank << " " << *distortion << std::endl;
}

void LayerVQ::CSL_ParallelCompetitiveLearning()
{
	//This is where we add the first touch of magic 
	//All PU:s has a copy of the training set thus 
	//Select a series of input vectors to perform CL on. 
	//n_s number of succesive inputs to each PU
	
	int n_s = m_x.size()/m_mpiSize;

	//effective memory management if placced un obect construct desruct. 
	//if(!m_winners)
//		m_winners = new int[n_s];

	vector<float> data;

	int max = (m_mpiRank*n_s)+n_s;
	if(m_mpiRank == m_mpiSize-1)
		max = m_x.size();
	
	int nrPart = max - m_mpiRank*n_s;
	m_winners = vector<int>(nrPart);

	//Getting Winners
	for(int index = m_mpiRank*n_s; index < max; index++)
	{
		m_winners[index - m_mpiRank*n_s] = CSL_GetClosestCodeVectorIndex( m_x[index],m_c);
	}

	//BroadCasting Winning vectors winning Vectors 

	//if(!m_total_winners)	
	//	m_total_winners = new int[m_maxUsedNodes*n_s];
	m_total_winners = vector<int>(m_mpiSize*n_s);//m_x.size());//m_maxUsedNodes*n_s);

	MPI_Allgather(&m_winners[0],n_s,MPI_INT,&m_total_winners[0],n_s,MPI_INT,NETWORK_COMM_WORLD);

	// distribute the extras if not evenly divisible
	if(m_mpiRank == m_mpiSize-1)
	{
		int nrRem = nrPart - n_s;
		MPI_Bcast(&nrRem,1,MPI_INT,m_mpiSize-1,NETWORK_COMM_WORLD);

		vector<int> extras(nrRem);

		for(int i=0;i<extras.size();i++)
		{
			extras[i] = m_winners[n_s+i];
			m_total_winners.push_back(extras[i]);
		}
		
		if(nrRem>0)
		{
			MPI_Bcast(&extras[0],nrRem,MPI_INT,m_mpiSize-1,NETWORK_COMM_WORLD);
		}
	}
	else
	{
		int nrRem;
		MPI_Bcast(&nrRem,1,MPI_INT,m_mpiSize-1,NETWORK_COMM_WORLD);
		if(nrRem>0)
		{
			vector<int> extras(nrRem);
			MPI_Bcast(&extras[0],nrRem,MPI_INT,m_mpiSize-1,NETWORK_COMM_WORLD);

			for(int i=0;i<extras.size();i++)
				m_total_winners.push_back(extras[i]);
		}
	}

	vector<int> order(m_x.size());

	if(m_mpiRank == 0)
	{
		for(int i=0;i<order.size();i++)
			order[i] = i;

		random_shuffle(order.begin(),order.end());
	}

	MPI_Bcast(&order[0],m_x.size(),MPI_INT,0,NETWORK_COMM_WORLD);

	//Adjust Code Book on accourding to vector of winners
	//if(!m_mpiRank)
	//	std::cout << "Adjusting Code Book According to winning indicies" << std::endl;
	if(m_c.size()<1)
	{
		cout<<m_mpiRank<<": m_c.size<1\n";
		cout.flush();
	}

	for (int i = 0; i < m_x.size();i++)//m_maxUsedNodes*n_s; i++ )
	{
		m_c[ m_total_winners[order[i]] ] = CSL_NewCodeVector( m_c[ m_total_winners[order[i]] ], m_x[ order[i] ] );
	}      	
}

int LayerVQ::CSL_GetClosestCodeVectorIndex(const vector<float>& data, const vector<vector<float> >& codeVectors)
{
  float minValue = 1e10;
  int minIndex = -1;
  
  for(int i=0;i<codeVectors.size();i++) {
    float val = CSL_RRValue(data,codeVectors[i]);//m_x[m],m_c[n]);
    if(val<minValue) {
      minValue = val;
      minIndex = i;
    }
  }
  
  return minIndex;
}


///////////////////////////////////////////////////////////////////////////////
/////// Connectivity
///////////////////////////////////////////////////////////////////////////////


void ProjectionModifierVQ::Modify()
{
	TimingStart(m_name);

	m_layerVQ->ModifyProjections(m_projection);
	
	TimingStop(m_name);
}

ProjectionModifierVQ::ProjectionModifierVQ(LayerVQ* layerVQ)
{
	m_eventId = 5; // unique id
	m_layerVQ = layerVQ;
	m_name = "ProjectionModifierVQ";
}

ProjectionModifierVQ* LayerVQ::GetProjectionModifier()
{
	return m_eventProjectionVQ;
}


// replace AddProjection to AddProjections etc
void LayerVQ::ModifyProjections(Projection* Projections)
{
	if(IsOn())
	{
		if(m_type == LayerVQ::VQStandard)
			m_totalWinnersUnits = m_population->MPI()->MPILayerGatherVariables(m_winnersUnits);

		PopulationColumns* post = (PopulationColumns*)Projections->PostLayer();
		PopulationColumns* pre = (PopulationColumns*)Projections->PreLayer();

		// rebuild connectivity
		Projections->Clear();
		vector<Unit*> postUnits = post->GetLocalRateUnits();

		vector<int> attachedHypercolumns;
		
		for(int i=0;i<postUnits.size();i++)
		{
			vector<long> preIds;
			vector<Hypercolumn*> postHc = ((RateUnit*)postUnits[i])->GetHypercolumns();

			// - this index should be matched to the index in totalWinnerUnits
			int hypercolumnIndex = postHc[0]->GetUnitIdLocal(); // assuming mc belong to only one hc

			attachedHypercolumns = m_totalWinnersUnits[hypercolumnIndex];
			if(attachedHypercolumns.size()>0)
			{
				vector<Hypercolumn*> hypercolumns = pre->GetHypercolumns();

				// check if Projections are between hypercolumns or minicolumns
				// fix
				
				int nrWinners = 0;
				//if(m_hasCheckedPreProjections == false)
				//{
					
					for(int i=0;i<m_totalWinnersUnits.size();i++)
						nrWinners+=m_totalWinnersUnits[i].size();

					//if(hypercolumns.size() != nrWinners) // always true atm
						m_isRateUnitPreProjections = true;

					m_hasCheckedPreProjections = true;
				//}

				vector<RateUnit*> minicolumns;
				vector<long> minicolumns2;

				if(m_isRateUnitPreProjections == true)
				{
					minicolumns = pre->GetRateUnits();
					
					for(int j=0;j<attachedHypercolumns.size();j++)
					{
						if(nrWinners == pre->GetNrHypercolumns()) // built up from hypercolumns (or minicolumn in each hypercolumn)
						{
							vector<long> mcIds = pre->GetRateUnitsIds(attachedHypercolumns[j]);

							for(int m=0;m<mcIds.size();m++)
								minicolumns2.push_back(mcIds[m]);
						}
						else	// built up from minicolumns
						{
							minicolumns2.push_back(attachedHypercolumns[j] + pre->GetUnitsStartId());
						}

						//Projections->AddProjection(minicolumns[attachedHypercolumns[j]]->GetUnitId(),postUnits[i]->GetUnitId(),false);
					}

					Projections->AddProjections(minicolumns2,postUnits[i]->GetUnitId(), postUnits[i]->GetUnitIdLocal());
					for(int j=0;j<minicolumns2.size();j++)
						Projections->GetConnectivity()->SetWeightValues(minicolumns2[j],postUnits[i]->GetUnitId());
					// not setting the weight or delay values in synapses hash - these will by default be 0
				}
				else
				{
					for(int j=0;j<attachedHypercolumns.size();j++)
					{
						minicolumns = hypercolumns[attachedHypercolumns[j]]->GetRateUnits();//k]->GetRateUnits();

						for(int m=0;m<minicolumns.size();m++)
						{
							minicolumns2.push_back(minicolumns[m]->GetUnitId());
							//Projections->AddProjection(minicolumns[m]->GetUnitId(),postUnits[i]->GetUnitId(),false);
						}
					}

					Projections->AddProjections(minicolumns2,postUnits[i]->GetUnitId(),postUnits[i]->GetUnitIdLocal());
					for(int j=0;j<minicolumns2.size();j++)
						Projections->GetConnectivity()->SetWeightValues(minicolumns2[j],postUnits[i]->GetUnitId());
					// not setting the weight or delay values in synapses hash - these will by default be 0
				}
			}
		}

		vector<int> localIndexes = post->GetLocalHypercolumnIndexes();
		for(int i=0;i<localIndexes.size();i++)
		{
			//m_winnersUnits
		}

		// not used anymore
		//Projections->GenerateHashTables(); // move to specific location
	}
}

/*
vector<float> time(1);
		time[0] = this->network()->GetCurrentTimeStep();

		vector<vector<float> > data;
		data.push_back(time);

		for(int i=0;i<m_Xi.size();i++)
			data.push_back(m_Xi[i]);

		return data;//m_Xi;*/

vector<vector<float> > LayerVQ::GetValuesToRecord()
{
	vector<vector<float> > f;

	if(IsRecording("stress"))
	{
		vector<vector<float> > data = m_recordedValues; // pop back instead with iterator

		for(int i=0;i<m_recordedValues.size();i++)
		{
			m_recordedValues[i].clear();
		}
		m_recordedValues.clear();

		return data;
	}
	else if(IsRecording())
	{
		for(int i=0;i<m_totalWinnersUnits.size();i++)
		{
			vector<float> f2;
			for(int j=0;j<m_totalWinnersUnits[i].size();j++)
			{
				f2.push_back((float)m_totalWinnersUnits[i][j]);
			}

			f.push_back(f2);
		}
	}

	return f;
}


///////////////////////////////////////////////////////////////////////////////
/////// CSL Implementation
///////////////////////////////////////////////////////////////////////////////

/*void CSL::Reset()
{
	m_M = 0;

}*/

void CSL::Initialize(vector<vector<float> >* x, int nrCodeVectors, int mpiRank, int mpiSize, bool printOutResult, NetworkObject* parent, MPI_Comm comm, int selectionNrs)
{
	m_parent = parent;
	m_useSelection = true;
	m_printOutResult = printOutResult;
	//m_x = *m_childPopulationModifier[0]->GetOutput(); // latest copy of (mds) data
	m_x = x;
	m_nrCodeVectors = nrCodeVectors;
	m_mpiRank = mpiRank;
	m_mpiSize = mpiSize;
	m_comm = comm;
	m_prevDistortion=0;
	m_distUnchanged=0;
	//m_epsilon = 0.01;//0.001;
	//m_eta = 0.01;//0.0001;//0.5;//0.001;
	
	// set initial values
	m_M = 1;
	m_currentDistortion = 1e10;
	//int mpiRank = this->m_population->network()->MPIGetNodeId();
	
	// randomly select nrNeurons different code vectors

	vector<int> selected;
	vector<int> indexes(m_x->size());

	for(int i=0;i<m_x->size();i++)
	{
		indexes[i] = i;
	}

	unsigned int seed;
	if(mpiRank == 0)
	{
		seed = (unsigned)time(NULL);
		if(m_comm != NULL)
			MPI_Bcast(&seed,1,MPI_UNSIGNED,0,m_comm); // should not be world
	}
	else
	{
		if(m_comm != NULL)
			MPI_Bcast(&seed,1,MPI_UNSIGNED,0,m_comm);
	}

	//srand(seed);//(unsigned)time(NULL));
	// randomization also in moving the codevectors (selection step)

	if(mpiRank == 0 && m_printOutResult)
	{
		cout<<mpiRank<<": m_x.size()="<<m_x->size()<<"\n";
		cout<<mpiRank<<": indexes.size()="<<indexes.size()<<"\n";
		cout<<mpiRank<<": m_nrCodeVectors="<<m_nrCodeVectors<<"\n";
		cout.flush();
	}

	for(int i=0;i<m_nrCodeVectors;i++)
	{
		// select indexes
		int r = 0;//rand() % indexes.size();
		if(indexes.size()==0)
			break;

		selected.insert(selected.end(),indexes[r]);
		indexes.erase(indexes.begin() + r);	
	}

 	// select as initial code vectors
	for(int i=0;i<selected.size();i++)
	{
		m_c.insert(m_c.end(),(*m_x)[selected[i]]);
	}

	if(mpiRank == 0 && m_printOutResult)
	{
		cout<<"CSL_Initiate complete.\n";
	}

	// Selection sizes (declining, determine number by input selectionNrs, default 3)
	for(int i=m_c.size()-m_c.size()/2; i>0; i--)
	{
		m_S.insert(m_S.end(),i);

		for(int m=0;m<selectionNrs;m++)
		{
			if(m_S.size()<m_x->size())
				m_S.insert(m_S.end(),i);
		}
		/*if(m_S.size()<m_x->size())
		m_S.insert(m_S.end(),i);
		if(m_S.size()<m_x->size())
		m_S.insert(m_S.end(),i);*/
	}

	m_M = 0; // Selection time step index
}

vector<float> CSL::NewCodeVector(vector<float> oldCodeVector, vector<float> inputVector)
{
	vector<float> newVector(oldCodeVector.size());

	//float eta = 0.01;
	for(int i=0;i<oldCodeVector.size();i++) {
		newVector[i] = oldCodeVector[i] + m_eta*(inputVector[i] - 
			oldCodeVector[i]);
	}

	return newVector;
}

// application specific: sensor drift
vector<float> CSL::NewCodeVectorVector(vector<float> oldCodeVector, vector<float> inputVector, vector<float> eta)
{
	vector<float> newVector(oldCodeVector.size());

	//float eta = 0.01;
	for(int i=0;i<oldCodeVector.size();i++) {
		newVector[i] = oldCodeVector[i] + eta[i]*(inputVector[i] - 
			oldCodeVector[i]);
	}

	return newVector;
}

void CSL::Selection(vector<int> indexes, vector<float> D)
{
	int s = indexes.size();

	// 1. Compute fitness measure using subdistortion
	float gamma = 0.5;

	vector<float> g(s);

	float sum = 0;
	for(int i=0;i<s;i++)
	{
		sum += pow(D[indexes[i]],gamma);
	}

	for(int i=0;i<s;i++)
	{
		g[i] = pow(D[indexes[i]],gamma) / sum;
	}

	// 2. Determine nr of neurons to be reproduced
	// i) compute integer part

	vector<int> u(s);
	for(int i=0;i<s;i++)
	{
		u[i] = (int)(g[i] * (float)s);
	}

	// ii) for top s - sum(u) wrt gs - u , add one

	int usum = 0;
	for(int i=0;i<u.size();i++) {
		usum += u[i];
	}

	vector<float> data(u.size());
	for(int i=0;i<u.size();i++) {
		data[i] = g[i]*(float)s - (float)u[i];
	}

	int nrValues = s - usum;

	vector<int> top = GetHighestIndexes(false, nrValues,data);

	// add one
	for(int i=0;i<top.size();i++) {
		u[top[i]] = u[top[i]] + 1;
	}

	// 3. For every j, reproduce u_j neurons by adding random perturbations (and remove if u_j == 0)
	vector<vector<float> > newC;

	/*
	unsigned int seed = GetSeed();

	srand(seed);//(unsigned)time(NULL));
	*/
	float epsilonMove = 0.01;

	for(int i=0;i<u.size();i++) 
	{
		//  m_c.erase(indexes[i]);

		for(int j=0;j<u[i];j++) 
		{
			vector<float> c(m_c[0].size());

			for(int x=0;x<m_c[0].size();x++) 
			{
				float dRand = ((float)rand()/(float)RAND_MAX);

				c[x] = (1+epsilonMove*(dRand - 0.5))*m_c[indexes[i]][x];
			}

			newC.insert(newC.end(),c);
		}
	}

	//TRACK: newC empty which raises exceptions from this Loop
	for(int i=0;i<indexes.size();i++) 
	{
		for(int j=0;j<m_c[0].size();j++)
			m_c[indexes[i]][j] = newC[i][j];
	}
}

vector<vector<float> > CSL::GetClosestCodeVectorMultipleIndexesAndValues(const vector<float>& data, const vector<vector<float> >& codeVectors, int maxSize)
{
  float minValue = 1e10;
  int minIndex = -1;

  vector<vector<float> > indexesAndValues;
  //vector<float> values;

  for(int i=0;i<codeVectors.size();i++) 
  {
    float val = RRValue(data,codeVectors[i]);//m_x[m],m_c[n]);

	int count = indexesAndValues.size();
	if(count == 0)
	{
		vector<float> iv(2);
		iv[0] = i;
		iv[1] = val;
		indexesAndValues.push_back(iv);
	}
	else 
	{
		bool placed = false;
		for(int j=0;j<count;j++)
		{
			if(val<indexesAndValues[j][1])
			{
				vector<float> iv(2);
				iv[0] = i;
				iv[1] = val;
				indexesAndValues.insert(indexesAndValues.begin()+j,iv);
				placed = true;
				break;
			}
		}

		if(placed == false)
		{
			if(count<maxSize+1)
			{
				vector<float> iv(2);
				iv[0] = i;
				iv[1] = val;
				indexesAndValues.push_back(iv);
			}
		}

		if(count>maxSize)
		{
			indexesAndValues.erase(indexesAndValues.begin()+maxSize,indexesAndValues.end());
		}
	}
  }
  
  return indexesAndValues;
}


vector<int> CSL::GetHighestIndexes(bool highestAndLowest, int nrValues, vector<float> data)
{
	vector<float> sorted;
	vector<int> indexesSorted;

	// sort
	for(int i=0;i<data.size();i++) 
	{
		bool placed = false;
		for(int j=0;j<sorted.size();j++) 
		{
			if(data[i]>sorted[j]) 
			{
				sorted.insert(sorted.begin() + j,data[i]);
				indexesSorted.insert(indexesSorted.begin() + j,i);
				placed = true;
				break;
			}
		}

		if(placed==false) {
			sorted.insert(sorted.end(),data[i]);
			indexesSorted.insert(indexesSorted.end(),i);
		}
	}

	// get nrValues first ones

	//Get the highest values
	if(highestAndLowest == false) {
		vector<int> out(nrValues);

		for(int i=0;i<nrValues;i++) {
			out[i] = indexesSorted[i];
		}

		return out;
	} else {
		// get alternatively highest and lowest values
		vector<int> out(nrValues);

		int indexHigh = 0;
		int indexLow = data.size()-1;//nrValues-1;

		for(int i=0;i<nrValues;i++) {
			if(i%2 == 0) { //even 
				out[i] = indexesSorted[indexHigh];
				indexHigh++;
			} else {
				out[i] = indexesSorted[indexLow];
				indexLow--;
			}
		}

		return out;
	}
}

unsigned int CSL::GetSeed()
{
	bool VQSeedUseTime = false;
	if(VQSeedUseTime == true)
	{
		int seed;

		if(m_mpiRank == 0)
		{
			seed = (unsigned)time(NULL);
			if(m_comm != NULL)
				MPI_Bcast(&seed,1,MPI_UNSIGNED,0,m_comm);
		}
		else
		{
			if(m_comm != NULL)
				MPI_Bcast(&seed,1,MPI_UNSIGNED,0,m_comm);
		}

		return seed;
	}
	else
	{
		return 0;
	}
}

void CSL::DriftAdaptStep(vector<float> x, float eta)
{
	// not parallel
	float oldEta = m_eta;
	m_eta = eta;//m_driftAdapt;

	int index = GetClosestCodeVectorIndex( x, m_c);
	m_c[ index ] = NewCodeVector( m_c[ index ], x );

	m_eta = oldEta;
}

void CSL::DriftAdaptStepVector(vector<float> x, vector<float> eta)
{
	int index = GetClosestCodeVectorIndex( x, m_c);
	m_c[ index ] = NewCodeVectorVector( m_c[ index ], x, eta );
}

vector<vector<int> > CSL::GetTotalWinnersUnits(int nrOverlaps)
{
	int n_s = m_x->size()/m_mpiSize;

	int max = (m_mpiRank*n_s)+n_s;
	if(m_mpiRank == m_mpiSize-1)
		max = m_x->size();
	
	int nrPart = max - m_mpiRank*n_s;
	vector<vector<int> > winners = vector<vector<int> >(nrOverlaps,vector<int>(nrPart));
	
	//Getting Winners
	for(int index = m_mpiRank*n_s; index < max; index++)
	{
		vector<vector<float> > indexesVals = GetClosestCodeVectorMultipleIndexesAndValues((*m_x)[index],m_c,nrOverlaps);
		for(int i=0;i<nrOverlaps;i++)
			winners[i][index - m_mpiRank*n_s] = (int)(indexesVals[i][0]);
		//m_winners[index - m_mpiRank*n_s] = GetClosestCodeVectorIndex( (*m_x)[index],m_c);
	}

	vector<vector<int> > total_winners = vector<vector<int> >(nrOverlaps,vector<int>(m_mpiSize*n_s));

	if(m_comm != NULL)
	{
		for(int i=0;i<nrOverlaps;i++)
		{
			MPI_Allgather(&winners[i][0],n_s,MPI_INT,&total_winners[i][0],n_s,MPI_INT,m_comm);
		}
	}
	else total_winners = winners;

	// distribute the extras if not evenly divisible
	if(m_mpiRank == m_mpiSize-1)
	{
		int nrRem = nrPart - n_s;
		if(m_comm != NULL)
			MPI_Bcast(&nrRem,1,MPI_INT,m_mpiSize-1,m_comm);

		vector<vector<int> > extras(nrOverlaps,vector<int>(nrRem));

		if(nrRem>0)
		{
			for(int j=0;j<nrOverlaps;j++)
			{
				for(int i=0;i<extras[j].size();i++)
				{
					extras[j][i] = winners[j][n_s+i];
					total_winners[j].push_back(extras[j][i]);
				}
			}


			if(m_comm != NULL)
			{
				for(int j=0;j<nrOverlaps;j++)
				{
					MPI_Bcast(&extras[j][0],nrRem,MPI_INT,m_mpiSize-1,m_comm);
				}
			}
		}
	}
	else
	{
		int nrRem;
		if(m_comm != NULL)
			MPI_Bcast(&nrRem,1,MPI_INT,m_mpiSize-1,m_comm);

		if(nrRem>0)
		{
			vector<vector<int> >extras(nrOverlaps,vector<int>(nrRem));
			if(m_comm != NULL)
			{
				for(int j=0;j<nrOverlaps;j++)
				{
					MPI_Bcast(&extras[j][0],nrRem,MPI_INT,m_mpiSize-1,m_comm);
				}
			}

			for(int i=0;i<extras.size();i++)
				for(int j=0;j<extras[i].size();j++)
					total_winners[i].push_back(extras[i][j]);
		}
	}

	m_totalWinnersUnits = vector<vector<int> >(m_c.size(),vector<int>(0));

	for(int i=0;i<total_winners.size();i++)
	{
		for(int j=0;j<total_winners[i].size();j++)
			m_totalWinnersUnits[total_winners[i][j]].push_back(j);
	}

	return m_totalWinnersUnits;
}

bool CSL::Step(vector<vector<float> >* x)
{
	m_x = x;
	//m_x = *m_childPopulationModifier[0]->GetOutput(); // latest copy of (mds) data

	float distortion = 0;
	vector<float> D(m_c.size());

	//adjust vectors using competitive learning

	ParallelCompetitiveLearning();
	ParallelDistortionCalculation( &distortion, &D );
	m_currentDistortion = distortion;
	if(fabs(m_prevDistortion-m_currentDistortion)<VQ_EPS) {
		m_distUnchanged++;
	}
	else {
		m_prevDistortion=m_currentDistortion;
		m_distUnchanged=0;
	}
	if(m_distUnchanged>=VQ_NUM)
		return false;
	//3. Selection
	if(fabs(distortion) < EPS) // usually never happens unless m_x.size() == m_c.size()
	{
		if(m_mpiRank == 0 && m_printOutResult)
		{
			cout<<"\ndistortion == 0!";
			cout.flush();
		}
	}
	else
	{
		if(m_S[m_M]>1 && m_useSelection == true)
		{
			//get the m_s[m_M] highest/lowest indexes wrt subdistortion

			vector<int> highestDIndexes = GetHighestIndexes(true,m_S[m_M],D);
			Selection(highestDIndexes, D);

			//if(count % m_settings.cslpCount == 0)
			//{
				m_M+=1;//2
				//if(m_M<0)
				//	m_M = 0;
				if(m_M>m_S.size()-1)
					m_M = m_S.size()-3;
			//}
		}
	}

	m_currentDistortion = distortion; // not needed to be global atm.
	if(m_mpiRank == 0 && DEBUG_LEVEL > 2 && m_printOutResult) 
	{
		cout << "VQ-dtot="<<distortion<< "[m_S[m_M]="<<m_S[m_M]<<"] ";
		cout.flush();
	}

	// for connectivity
	m_totalWinnersUnits.clear();
	m_totalWinnersUnits = vector<vector<int> >(m_c.size(),vector<int>(0));

	for(int i=0;i<m_total_winners.size();i++)
	{
		m_totalWinnersUnits[m_total_winners[i]].push_back(i);
	}
	return true;
}

float CSL::RRValue(const vector<float>& x1, const vector<float>& x2)
{
	float val = 0;
	for(int i=0;i<x1.size();i++) 
	{
		val+=(x1[i]-x2[i])*(x1[i]-x2[i]);//pow(x1[i]-x2[i],2);
	}
	return sqrt(val);//val;
}

// Currently only saves/loads the values of the code vectors
void CSL::SaveState()
{
	if(m_c.size()>0)
	{
		m_cSaved = vector<vector<float> >(m_c.size(),vector<float>(m_c[0].size()));

		for(int i=0;i<m_c.size();i++)
			for(int j=0;j<m_c[i].size();j++)
				m_cSaved[i][j] = m_c[i][j];
	}
}

void CSL::LoadState()
{
	if(m_cSaved.size()>0)
	{
		m_c = vector<vector<float> >(m_cSaved.size(),vector<float>(m_cSaved[0].size()));
		
		for(int i=0;i<m_cSaved.size();i++)
			for(int j=0;j<m_cSaved[i].size();j++)
				m_c[i][j] = m_cSaved[i][j];

//		m_cSaved.clear();
	}
}

vector<float> CSL::GetClosestCodeVectorIndexAndValue(const vector<float>& data, const vector<vector<float> >& codeVectors)
{
	vector<float> out;

	float minValue = 1e10;
	float minIndex = -1;

	for(int i=0;i<codeVectors.size();i++) 
	{
		float val = RRValue(data,codeVectors[i]);//m_x[m],m_c[n]);

		if(val<minValue) 
		{
			minValue = val;
			minIndex = i;
		}
	}

	out.insert(out.end(),minIndex);
	out.insert(out.end(),minValue);

	return out;
}

vector<int> CSL::GetClosestCodeVectorIndexes(const vector<vector<float> >& data)
{
	vector<int> out(data.size());

	for(int i=0;i<data.size();i++)
	{
		out[i] = GetClosestCodeVectorIndex(data[i],m_c);
	}

	return out;
}


void CSL::ParallelDistortionCalculation(float *distortion, vector<float> *D)
{
	int n_s = m_x->size()/m_mpiSize;

	vector<float> totD(m_c.size());
	vector<float> subD(m_c.size());

	float newDistortion = 0;

	int max = (m_mpiRank*n_s)+n_s;
	if(m_mpiRank == m_mpiSize-1)
		max = m_x->size();

	for(int i = m_mpiRank*n_s; i < max; i++)
	{
		vector<float> v =  GetClosestCodeVectorIndexAndValue((*m_x)[i], m_c);
		subD[(int)v[0]] += v[1];
		newDistortion += v[1];
	}

	if(m_comm != NULL)
	{
 		MPI_Allreduce(&newDistortion, distortion, 1, MPI_FLOAT, MPI_SUM, m_comm);
		MPI_Allreduce(&subD[0], &totD[0], m_c.size(), MPI_FLOAT, MPI_SUM, m_comm);
	}
	else
	{
		*distortion = newDistortion;
		totD = subD;
	}

	for(int i=0;i<totD.size();i++)
		totD[i] = totD[i]/(float)m_x->size();

	*D = totD;
	*distortion = *distortion / (float)m_x->size();

	if(m_mpiRank == 0 && m_printOutResult)
	{
		cout<<"dist="<<*distortion<<" ";cout.flush();
	}
	//MPI_Barrier(NETWORK_COMM_WORLD);
	//std::cout << "Rank " << m_mpiRank << " " << *distortion << std::endl;
}

void CSL::ParallelCompetitiveLearning()
{
	//This is where we add the first touch of magic 
	//All PU:s has a copy of the training set thus 
	//Select a series of input vectors to perform CL on. 
	//n_s number of succesive inputs to each PU
	
	int n_s = m_x->size()/m_mpiSize;

	//effective memory management if placced un obect construct desruct. 
	//if(!m_winners)
//		m_winners = new int[n_s];

	vector<float> data;

	int max = (m_mpiRank*n_s)+n_s;
	if(m_mpiRank == m_mpiSize-1)
		max = m_x->size();
	
	int nrPart = max - m_mpiRank*n_s;
	m_winners = vector<int>(nrPart);

	//Getting Winners
	for(int index = m_mpiRank*n_s; index < max; index++)
	{
		m_winners[index - m_mpiRank*n_s] = GetClosestCodeVectorIndex( (*m_x)[index],m_c);
	}

	//BroadCasting Winning vectors winning Vectors 

	//if(!m_total_winners)	
	//	m_total_winners = new int[m_maxUsedNodes*n_s];
	m_total_winners = vector<int>(m_mpiSize*n_s);//m_x.size());//m_maxUsedNodes*n_s);

	if(m_parent!=NULL)
		m_parent->TimingStart("CommunicationCSL");

	if(m_comm != NULL)
		MPI_Allgather(&m_winners[0],n_s,MPI_INT,&m_total_winners[0],n_s,MPI_INT,m_comm);
	else m_total_winners = m_winners;

	// distribute the extras if not evenly divisible
	if(m_mpiRank == m_mpiSize-1)
	{
		int nrRem = nrPart - n_s;
		if(m_comm != NULL)
			MPI_Bcast(&nrRem,1,MPI_INT,m_mpiSize-1,m_comm);

		vector<int> extras(nrRem);

		for(int i=0;i<extras.size();i++)
		{
			extras[i] = m_winners[n_s+i];
			m_total_winners.push_back(extras[i]);
		}
		
		if(nrRem>0)
		{
			if(m_comm != NULL)
				MPI_Bcast(&extras[0],nrRem,MPI_INT,m_mpiSize-1,m_comm);
		}
	}
	else
	{
		int nrRem;
		if(m_comm != NULL)
			MPI_Bcast(&nrRem,1,MPI_INT,m_mpiSize-1,m_comm);

		if(nrRem>0)
		{
			vector<int> extras(nrRem);
			if(m_comm != NULL)
				MPI_Bcast(&extras[0],nrRem,MPI_INT,m_mpiSize-1,m_comm);

			for(int i=0;i<extras.size();i++)
				m_total_winners.push_back(extras[i]);
		}
	}

	vector<int> order(m_x->size());

	if(m_mpiRank == 0)
	{
		for(int i=0;i<order.size();i++)
			order[i] = i;

		random_shuffle(order.begin(),order.end());
	}

	if(m_comm != NULL)
		MPI_Bcast(&order[0],m_x->size(),MPI_INT,0,m_comm);

	if(m_parent!=NULL)
		m_parent->TimingStop("CommunicationCSL");

	//Adjust Code Book on accourding to vector of winners
	//if(!m_mpiRank)
	//	std::cout << "Adjusting Code Book According to winning indicies" << std::endl;
	if(m_c.size()<1)
	{
		cerr<<m_mpiRank<<": m_c.size<1\n";
		cout.flush();
	}

	// optimize (!)
	for (int i = 0; i < m_x->size();i++)//m_maxUsedNodes*n_s; i++ )
	{
		if(m_total_winners[order[i]] > m_c.size()-1)
		{
			cerr<<"Error tot winners size. ("<<m_total_winners[order[i]]<<", "<<m_c.size()<<")";cout.flush();
		}

		m_c[ m_total_winners[order[i]] ] = NewCodeVector( m_c[ m_total_winners[order[i]] ], (*m_x)[ order[i] ] );
	}      	
}

int CSL::GetClosestCodeVectorIndex(const vector<float>& data, const vector<vector<float> >& codeVectors)
{
	float minValue = 1e10;
	int minIndex = -1;
	float val;

	for(int i=0;i<codeVectors.size();i++) 
	{
		val = RRValue(data,codeVectors[i]);//m_x[m],m_c[n]);
		if(val<minValue) 
		{
			minValue = val;
			minIndex = i;
		}
	}

	if(minIndex == -1)
	{
		cout<<minValue<<"/"<<val<<"";
		for(int i=0;i<data.size();i++)
		{
			cout<<data[i]<<" ";
		}
		cout<<";";
		cout.flush();
	}

	return minIndex;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// ProjectionModifierCSL: 
//////////////// CSL usage as feature extractor in multiple hypercolumns (simultaneously, in parallel)
///////////////////////////////////////////////////////////////////////////////////////////////////////////


ProjectionModifierCSL::ProjectionModifierCSL()
{
	m_eventId = 9;

	m_eta = 0.01f;
	m_epsilon = 0.01;

	m_transferFunction = new TransferCSL();
	m_maxPatterns = 0;
	m_name = "ProjectionModifierCSL";
}

void ProjectionModifierCSL::Initialize(Projection* Projection)
{
	network(Projection->network());
	vector<float> postValues = m_projectionFixed->GetPostValues();

	m_projectionFixed = Projection;
	m_idsPost = m_projectionFixed->GetPostIds();

	for(int i=0;i<m_idsPost.size();i++)
	{
		int localIndex = ((RateUnit*)network()->GetUnitFromId(m_idsPost[i]))->GetUnitIdLocalInHypercolumn();
		m_idsIndexInHc[m_idsPost[i]] = localIndex;
	}

	// build communicators for output hypercolumns
	m_projectionFixed->PostLayer()->MPI()->MPICreateCommsHypercolumns();
	vector<MPI_Comm*> commsAll = m_projectionFixed->PostLayer()->MPI()->MPIGetCommHCs(); // should locally only create the needed comms

	bool allSame = true;
	for(int i=1;i<commsAll.size();i++)
	{
		if((commsAll[i]) != commsAll[0])
		{
			allSame = false;
			break;
			//cout<<"(E) commAll["<<i<<"]==NULL";cout.flush();
		}
	}

	if(allSame == true)
	{
		cout<<"(E) commsAll the same.("<<commsAll[0]<<")";
	}

	vector<int> hcIndexes; // indexes of the output hcs
	vector<int> hcAllIndexes; // index of each mc in hc - which it belongs to (in m_csl)
	vector<int> nrCodeVectors;
	vector<long> unitIdReprHc;

	for(int i=0;i<m_idsPost.size();i++)
	{
		int hcId = ((RateUnit*)network()->GetUnitFromId(m_idsPost[i]))->GetHypercolumnId();
		int hcIndex = ((Hypercolumn*)network()->GetUnitFromId(hcId))->GetUnitIdLocal();//m_projectionFixed->PostLayer()->

		bool contains = false;
		int index=-1;
		for(int j=0;j<hcIndexes.size();j++)
		{
			if(hcIndexes[j] == hcIndex)
			{
				contains = true;
				index = j;
				break;
			}
		}

		if(contains == false)
		{
			unitIdReprHc.push_back(m_idsPost[i]);
			hcIndexes.push_back(hcIndex);
			if(commsAll.size()-m_comms.size()>0)
				m_comms.push_back(commsAll[((Hypercolumn*)network()->GetUnitFromId(hcId))->GetUnitIdLocal()]);//m_comms.size()]);// commsAll[hcIndexes.size()-1]); //this ok?
				//m_comms.push_back(commsAll[hcIndex]);
			hcAllIndexes.push_back(1);//hcIndex);
			//nrCodeVectors.push_back(m_projectionFixed->PostLayer()->GetStructure()[hcIndex]);
			nrCodeVectors.push_back(((Hypercolumn*)network()->GetUnitFromId(hcId))->GetTotalNrRateUnits()); // nr minicolumns in hypercolumn
		}
		else
			hcAllIndexes[hcAllIndexes.size()-1] = hcAllIndexes[hcAllIndexes.size()-1] + 1;//hcIndexes[index]);// could be 'index' instead
	}

	for(int i=0;i<hcIndexes.size();i++)
	{
		CSL* csl = new CSL(m_epsilon,m_eta);
		csl->SetMPIParameters(this->network()->MPIGetNodeId(), this->network()->MPIGetNrProcs());
		m_csl.push_back(csl);
		vector<vector<float> > data;
		m_x.push_back(data);

		vector<long> processIds = ((PopulationColumns*)m_projectionFixed->PostLayer())->GetNodeIndexes(hcIndexes[i]);
		//vector<int> locSize(processIds.size());
		//for(int m=0;m<locSize.size();m++) locSize[m] = processIds[m];
		m_mpiLocSize.push_back(processIds.size());
		int mpiRank = this->network()->MPIGetNodeId();
//		cout<<processIds.size()<<" ";cout.flush();
		m_mpiLocRank.push_back(mpiRank-processIds[0]); // assume a linear distribution of nodes
		m_firstRun.push_back(true);
	}

	m_mcReprHc = unitIdReprHc;
	m_nrCodeVectors = nrCodeVectors;
	m_hcIndexes = hcIndexes;
	m_hcAllIndexes = hcAllIndexes;
}

void ProjectionModifierCSL::SetProjection(Projection* c)
{
	m_projectionFixed = (ProjectionFixed*)c;
}

void ProjectionModifierCSL::SaveState()
{
	for(int i=0;i<m_csl.size();i++)
	{
		m_csl[i]->SaveState();
	}
}

void ProjectionModifierCSL::LoadState()
{
	for(int i=0;i<m_csl.size();i++)
	{
		m_csl[i]->LoadState();
	}
}

void ProjectionModifierCSL::DriftAdaptStep(float eta)
{
	// not necessary to put in m_x in current impl.
	for(int i=0;i<m_csl.size();i++)
	{
		if(m_maxPatterns!=0)
		{
			if(m_x[i].size()>=m_maxPatterns)
			{
				m_x[i].erase(m_x[i].begin());// remove earliest added
			}

			vector<float> preValues = m_projectionFixed->GetPreValues(m_mcReprHc[i]);

			if(preValues.size()>0)
			{
				m_x[i].push_back(preValues);
			}

			m_csl[i]->DriftAdaptStep(preValues, eta);
		}
	}
}

void ProjectionModifierCSL::DriftAdaptStepVector(vector<float> eta)
{
	// not necessary to put in m_x in current impl.
	for(int i=0;i<m_csl.size();i++)
	{
		if(m_maxPatterns!=0)
		{
			if(m_x[i].size()>=m_maxPatterns)
			{
				m_x[i].erase(m_x[i].begin());// remove earliest added
			}

			vector<float> preValues = m_projectionFixed->GetPreValues(m_mcReprHc[i]);

			if(preValues.size()>0)
			{
				m_x[i].push_back(preValues);
			}

			m_csl[i]->DriftAdaptStepVector(preValues, eta);
		}
	}
}

void ProjectionModifierCSL::Modify()
{
	TimingStart(m_name);

	if(IsOn() == false) return;

	int processId = m_projectionFixed->PreLayer()->network()->MPIGetNodeId(); // will be put in mpi-class

	if(processId == 0)
	{
		cout<<".";
		cout.flush();
	}

	//vector<float> postValues = m_projectionFixed->GetPostValues();
	
	vector<vector<long> >* preIds;// = m_projectionFixed->PreIds(); // move to initializer (until Invalidate().. called)

	long preId, postId;
	float weight;
	double totWeights;
	float dw;


	int currentIndex = 0;
	for(int i=0;i<m_csl.size();i++) // go through all csl:s this process is involved with
	{
		if(m_maxPatterns!=0)
		{
			if(m_x[i].size()>=m_maxPatterns)
			{
				m_x[i].erase(m_x[i].begin());// remove earliest added
			}

			vector<float> preValues = m_projectionFixed->GetPreValues(m_mcReprHc[i]);

			if(preValues.size()>0)
			{
				m_x[i].push_back(preValues);

				bool hasRun = false;
				if(m_firstRun[i] == true)
				{
					// set up a flag to this (will be as a diff between iterative and batch)
					if(m_x[i].size() >= m_maxPatterns)
					//if(m_x[i].size() >= m_nrCodeVectors[i])
					{
						bool printOut = false;
						if(this->network()->MPIGetNodeId() == 0)
							printOut = true;

						if(m_comms[i] == NULL)
							m_csl[i]->Initialize(&(m_x[i]),m_nrCodeVectors[i],m_mpiLocRank[i],m_mpiLocSize[i],printOut,this,NULL);
						else
							m_csl[i]->Initialize(&(m_x[i]),m_nrCodeVectors[i],m_mpiLocRank[i],m_mpiLocSize[i],printOut,this, *(m_comms[i]));//this->network()->MPIGetNodeId(),this->network()->MPIGetNrProcs(),*(m_comms[i]));
						
						m_csl[i]->Step(&(m_x[i]));
						m_firstRun[i] = false;
						hasRun = true;
					}
				}
				else
				{
					if(m_x[i].size()!=m_csl[i]->GetCodeVectors()->size()) // could put check inside csl
						SwitchOnOff(m_csl[i]->Step(&(m_x[i])));
					
					hasRun = true;
				}

				if(hasRun == true)
				{
					// set the incoming weights to the code vectors
					vector<vector<float> >* cvs = m_csl[i]->GetCodeVectors();

					currentIndex = 0;
					for(int m=0;m<i;m++)
						currentIndex+=m_hcAllIndexes[m];

					for(int j=0;j<m_hcAllIndexes[i];j++)
					{
						long postId = m_idsPost[currentIndex];
						long preId;

						int indexInHc = m_idsIndexInHc[postId]; // to determine which codevector should be assigned to the unit

						vector<long> preIdsPart =  m_projectionFixed->GetPreIds(postId);//m_projectionFixed->GetPreIds(postId); // slower method

						if(preIdsPart.size()>0) // ! ?
							for(int k=0;k<preValues.size();k++)//preIdsPart.size();k++)
							{
								//preId = (*preIds)[currentIndex][k];
								preId = preIdsPart[k];
								network()->SetWeight((*cvs)[indexInHc][k],preId,postId);
							}

							currentIndex++;
					}
				}
			}
		}
	}
	
	// add pre values to m_x (or replace if too large)
	// run 1 new step

	/*for(int j=0;j<postValues.size();j++)
	{	
		if( postValues[j] == 1.0) // assume WTA(1.0) is used
		{
			postId = m_idsPost[j];

			vector<float> preValues = m_projectionFixed->GetPreValues(m_idsPost[j]);

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

		((RateUnit*)m_projectionFixed->PostLayer()->network()->GetUnitFromId(m_idsPost[j]))->AddInhibBeta(m_b[j]);
	}*/

	TimingStop(m_name);
}

// assumes only one csl instance?
void ProjectionModifierCSL::UpdateWeights(int i, int currentIndex)
{
	// set the incoming weights to the code vectors
	vector<vector<float> >* cvs = m_csl[i]->GetCodeVectors();
	vector<float> preValues = m_projectionFixed->GetPreValues(m_mcReprHc[i]);

	for(int j=0;j<m_hcAllIndexes[i];j++)
	{
		long postId = m_idsPost[currentIndex];
		long preId;

		int indexInHc = m_idsIndexInHc[postId]; // to determine which codevector should be assigned to the unit

		vector<long> preIdsPart =  m_projectionFixed->GetPreIds(postId);//m_projectionFixed->GetPreIds(postId); // slower method

		if(preIdsPart.size()>0) // ! ?
			for(int k=0;k<preValues.size();k++)
			{
				//preId = (*preIds)[currentIndex][k];
				preId = preIdsPart[k];
				network()->SetWeight((*cvs)[indexInHc][k],preId,postId);
			}

			currentIndex++;
	}
}

void ProjectionModifierCSL::Simulate(UnitModifier* e)
{
}