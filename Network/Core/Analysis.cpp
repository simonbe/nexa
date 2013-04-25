#include <cmath>
#include "Network.h"
#include "Analysis.h"
#include "NetworkPopulation.h"
#include "MPIDistribution.h"

void Analysis::DoOwnFileWriting()
{
	if(m_useOwnFileWriting == true && this->m_network->MPIGetNodeId()==0)
	{
		if(m_results.size()>0)
		{
			cout<<"Own File Writing...";cout.flush();
			/////////////////////////////////////////////
			// do own file writing if (large) matrix should not be in same analysis file as others
			////////////////////////////////////////////

			Storage storage;
			storage.SaveDataFloatCSV((char*)(m_ownFilename.c_str()),m_results,Storage::Standard,"");

			// return empty to avoid writing anything more to disk
			/*for(int i=0;i<m_results.size();i++)
				m_results[i].clear();

			m_results.clear();*/
			cout<<"ok.";cout.flush();
		}
	}
	else if(this->m_network->MPIGetNodeId()!=0)
		m_results.clear();
	//	else m_results.clear();
}

AnalysisLayer::AnalysisLayer(Population* Population) // use
{
	m_populations.push_back(Population);
	m_name = "AnalysisLayer";
	this->m_network = Population->network();
}

AnalysisLayer::AnalysisLayer(vector<Population*> Populations) // use
{
	m_populations = Populations;
	m_name = "AnalysisLayer";
	m_on = true;
	this->m_network = Populations[0]->network();
}

FisherRatio::FisherRatio(Population* Population)
{
	m_populations.push_back(Population);
	m_useIncrementalVersion = false;
	m_name = "FisherRatio";
	this->m_network = Population->network();
}

FisherRatio::FisherRatio(Population* Population, bool useIncrementalVersion)
{
	m_populations.push_back(Population);
	m_populations.push_back(Population);
	m_useIncrementalVersion = useIncrementalVersion;
	m_name = "FisherRatio";
	this->m_network = Population->network();
}

FisherRatio::FisherRatio(Population* Population1, Population* Population2, bool useIncrementalVersion)
{
	m_populations.push_back(Population1);
	m_populations.push_back(Population2);
	m_useIncrementalVersion = useIncrementalVersion;
	m_name = "FisherRatio";
	this->m_network = Population1->network();
}

void AnalysisTiming::Finalize()
{
	// calculate mean across processes

	if(m_timingTot.size()>0)
	{
		// find out the maximum number
		int max = -1;
		int thisNr = m_timingTot.size();
		MPI_Allreduce(&thisNr,&max,1,MPI_INT,MPI_MAX,NETWORK_COMM_WORLD);

		map<string,double>::iterator p;
		vector<float> totValues(max);
		vector<float> sqrmValues(max);
		vector<float> sqrValues(max);
		m_meanValues = vector<float>(max);
		m_stdValues = vector<float>(max);
		m_locMin = vector<int>(max);
		m_locMax = vector<int>(max);
		m_min = vector<float>(max);
		m_max = vector<float>(max);
				
		vector<float> maxValues(max);
		vector<float> maxSqrmValues(max);
		vector<float> maxSqrValues(max);
		m_meanMaxValues = vector<float>(max);
		m_stdMaxValues = vector<float>(max);
		vector<float> minValues(max);
		vector<float> minSqrmValues(max);
		vector<float> minSqrValues(max);
		m_meanMinValues = vector<float>(max);
		m_stdMinValues = vector<float>(max);

		int index = 0;
		for(p = m_timingTot.begin(); p != m_timingTot.end(); p++)
		{
			totValues[index] = p->second;
			index++;
		}
			
		index = 0;
		for(p = m_timingMax.begin(); p != m_timingMax.end(); p++)
		{
			maxValues[index] = p->second;
			index++;
		}
		
		index = 0;
		for(p = m_timingMin.begin(); p != m_timingMin.end(); p++)
		{
			minValues[index] = p->second;
			index++;
		}

		if(thisNr!=max)
		{
			totValues = vector<float>(max);
			maxValues = vector<float>(max);
			minValues = vector<float>(max);
		}

		MPI_Allreduce(&totValues[0],&m_meanValues[0],max,MPI_FLOAT,MPI_SUM,NETWORK_COMM_WORLD);

		for(int i=0;i<m_meanValues.size();i++)
		{
			m_meanValues[i] = m_meanValues[i]/(float)m_network->MPIGetNrProcs();
		}

		for(int i=0;i<m_meanValues.size();i++)
		{
			sqrmValues[i] = (totValues[i]-m_meanValues[i])*(totValues[i]-m_meanValues[i]);//pow((totValues[i]-m_meanValues[i]),2);
		}

		// could do only reduce
		MPI_Allreduce(&sqrmValues[0],&sqrValues[0],max,MPI_FLOAT,MPI_SUM,NETWORK_COMM_WORLD);

		// std
		for(int i=0;i<m_meanValues.size();i++)
		{
			m_stdValues[i] = sqrt(sqrValues[i]/(float)m_network->MPIGetNrProcs());
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// same for the max values - will show variations in the timings (may lead to waiting times at barriers)

		MPI_Allreduce(&maxValues[0],&m_meanMaxValues[0],max,MPI_FLOAT,MPI_SUM,NETWORK_COMM_WORLD);
		MPI_Allreduce(&minValues[0],&m_meanMinValues[0],max,MPI_FLOAT,MPI_SUM,NETWORK_COMM_WORLD);

		for(int i=0;i<m_meanMaxValues.size();i++)
		{
			m_meanMaxValues[i] = m_meanMaxValues[i]/(float)m_network->MPIGetNrProcs();
			m_meanMinValues[i] = m_meanMinValues[i]/(float)m_network->MPIGetNrProcs();
		}

		for(int i=0;i<m_meanMaxValues.size();i++)
		{
			maxSqrmValues[i] = (maxValues[i]-m_meanMaxValues[i])*(maxValues[i]-m_meanMaxValues[i]);//pow((maxValues[i]-m_meanMaxValues[i]),2);
			minSqrmValues[i] = (minValues[i]-m_meanMinValues[i])*(minValues[i]-m_meanMinValues[i]);//pow((minValues[i]-m_meanMinValues[i]),2);
		}

		// could do only reduce
		MPI_Allreduce(&maxSqrmValues[0],&maxSqrValues[0],max,MPI_FLOAT,MPI_SUM,NETWORK_COMM_WORLD);
		MPI_Allreduce(&minSqrmValues[0],&minSqrValues[0],max,MPI_FLOAT,MPI_SUM,NETWORK_COMM_WORLD);

		// std
		for(int i=0;i<m_meanMaxValues.size();i++)
		{
			m_stdMaxValues[i] = sqrt(maxSqrValues[i]/(float)m_network->MPIGetNrProcs());
			m_stdMinValues[i] = sqrt(minSqrValues[i]/(float)m_network->MPIGetNrProcs());
		}

		// for minloc, maxloc
		struct {
			float value;
			int   index;
		} m_in, m_out;

		// min, max locations of mean values
		index = 0;
		for(int m=0;m<totValues.size();m++)//p = m_timingTot.begin(); p != m_timingTot.end(); p++)
		{
			m_in.index = m_network->MPIGetNodeId();
			m_in.value = totValues[m];//p->second;

			MPI_Allreduce(&m_in,&m_out,1,MPI_FLOAT_INT,MPI_MINLOC,NETWORK_COMM_WORLD);
			m_locMin[index] = m_out.index;
			m_min[index] = m_out.value;

			MPI_Allreduce(&m_in,&m_out,1,MPI_FLOAT_INT,MPI_MAXLOC,NETWORK_COMM_WORLD);
			m_locMax[index] = m_out.index;
			m_max[index] = m_out.value;

			index++;
		}						
	}
}

void FisherRatio::Simulate()
{
	if(!IsOn()) return;

	AnalysisLayer::Simulate();

	if(m_continue == true)
	{

		if(m_useIncrementalVersion == false) // batch version
		{ 
			// only add current pattern to current class
			vector<Unit*> localUnits;
			if(m_populations.size() == 1)
				localUnits = ((PopulationColumns*)m_populations[0])->GetUnits();
			else
				localUnits = ((PopulationColumns*)m_populations[this->m_currentPatternClass])->GetUnits();

			vector<float> values(localUnits.size());

			for(int i=0;i<localUnits.size();i++)
			{
				values[i] = localUnits[i]->GetValue();
			}

			for(int i=0;i<localUnits.size();i++)
			{
				values[i] = localUnits[i]->GetValue();
			}

			if(this->m_data.size() < this->m_currentPatternClass+1)
			{
				for(int i=m_data.size();i<this->m_currentPatternClass+1;i++)
				{
					m_data.push_back(vector<vector<float> >(0));
				}
			}

			this->m_data[this->m_currentPatternClass].push_back(values);
		}
		else
		{
			// incremental version, do as in NetworkCorr
		}
	}
}

void FisherRatio::CalcFisherRatio()
{
	if(m_useIncrementalVersion == false) // batch version
	{ 
		// calculate results

		// means
		vector<vector<float> > means(this->m_data.size());
		for(int i=0;i<m_data.size();i++)
		{
			means[i] = vector<float>(m_data[i].size());

			if(m_data.size()>0)
			{
				for(int j=0;j<m_data[i].size();j++)
				{
					for(int k=0;k<m_data[i][j].size();k++)
						means[i][j] += m_data[i][j][k];

					means[i][j] /= (float)m_data[i].size();
				}
			}
		}

		// dist sqr
		vector<vector<float> > localDist(this->m_data.size());
		for(int i=0;i<m_data.size();i++)
		{
			if(m_data[i].size()>0)
			{
				localDist[i] = vector<float>(m_data[i].size());

				for(int j=0;j<m_data[i].size();j++)
				{
					for(int k=0;k<m_data[i][j].size();k++)
						localDist[i][j] += (m_data[i][j][k]-means[i][j])*(m_data[i][j][k]-means[i][j]);//pow(m_data[i][j][k]-means[i][j],2);
				}
			}
		}

		// put together, sqrt
		vector<vector<float> > dist(this->m_data.size());
		for(int i=0;i<localDist.size();i++)
		{
			vector<float> f;
			if(localDist[i].size()>0)
				f = m_populations[0]->MPI()->MPILayerReduceVariable(localDist[i]); // assumes same mpi distribution for both populations (will be by default)

			dist[i] = f;
		}

		for(int i=0;i<dist.size();i++)
		{
			for(int j=0;j<dist[i].size();j++)
				dist[i][j]=sqrt(dist[i][j]);
		}

		// vars
		vector<float> vars(this->m_data.size());
		for(int i=0;i<vars.size();i++)
		{
			for(int j=0;j<dist[i].size();j++)
				vars[i] += dist[i][j];

			vars[i]/=(float)dist.size();
		}

		// fisher ratio between all classes

		m_results = vector<vector<float> >(means.size());
		for(int i=0;i<means.size();i++)
		{
			vector<float> ratios;
			for(int j=0;j<=i;j++)
			{
				float meanDist = 0;
				int minSize = min(means[i].size(),means[j].size()); // in case one is of zero size
				for(int m=0;m<minSize;m++)
				{
					meanDist+=(means[i][m]-means[j][m])*(means[i][m]-means[j][m]);//pow(means[i][m]-means[j][m],2);
				}

				meanDist = sqrt(meanDist);
				float ratio = -1;
				if(minSize>0)
					ratio = meanDist/(vars[i]+vars[j]);

				ratios.push_back(ratio);
			}

			m_results[i] = ratios;
		}
	}
	else // incremental version
	{
		// do nothing
	}
}

void FisherRatio::Finalize()
{
	CalcFisherRatio();
	DoOwnFileWriting();
}


void AnalysisLayer::Simulate()
{
	if(!IsOn()) return;

	m_continue = true;

	if(m_samplingTimeSteps>0)
	{
		m_currentSamplingTimeStep++;

		if(m_currentSamplingTimeStep<m_samplingStartTimeStep)
			m_continue = false;
		else
		{
			if((m_currentSamplingTimeStep-m_samplingStartTimeStep) % m_samplingTimeSteps == 0)
			{
			}
			else m_continue = false;
		}
	}

}

void AnalysisDistance::Simulate()
{
	AnalysisLayer::Simulate();
		
	if(m_continue ==  true)
	{
		vector<RateUnit*> localUnits = ((PopulationColumns*)m_populations[0])->GetRateUnits();//((PopulationColumns*)m_populations[0])->GetLocalRateUnits();
		vector<float> values(localUnits.size());

		for(int i=0;i<localUnits.size();i++)
		{
			values[i] = localUnits[i]->GetValue();
		}

		// collect this representation to test against at a later time
		if(this->m_collectRepresentations == true)
		{
			this->AddSpecificRepresentation(values);
		}
		else // compare representation to either previous timesteps or collected/set representations
		{
			bool doCalc = true;

			if(this->m_nrConcatenatedItems>0)
			{
				for(int i=0;i<values.size();i++)
					m_buffer2.push_back(values[i]);

				if(m_buffer2.size() == m_nrConcatenatedItems*values.size())
				{
					values = m_buffer2;
					m_buffer2.clear();
				}
				else
					doCalc = false;
			}

			if(doCalc == true)
			{
				if(m_specificRepresentations.size() == 0)
					m_buffer.push_back(values);
				else
					m_buffer = m_specificRepresentations; // higher memory usage than necessary, could get high if large training set (!)

				if(this->m_nrPreviousTimeSteps>0)
				{
					if(m_buffer.size()>m_nrPreviousTimeSteps)
						m_buffer.erase(m_buffer.begin(),m_buffer.begin()+1);
				}

				vector<float> distances(m_buffer.size());

				if(m_measure == AnalysisDistance::Abs)
				{
					for(int i=0;i<m_buffer.size();i++)
						for(int j=0;j<m_buffer[i].size();j++)
							distances[i]+=abs(m_buffer[i][j]-values[j]);
				}
				else if(m_measure == AnalysisDistance::Euclidean)
				{
					for(int i=0;i<m_buffer.size();i++)
						for(int j=0;j<m_buffer[i].size();j++)
							distances[i]+=(m_buffer[i][j]-values[j])*(m_buffer[i][j]-values[j]);//pow(m_buffer[i][j]-values[j],2);

				}
				else if(m_measure == AnalysisDistance::Hamming)
				{
					for(int i=0;i<m_buffer.size();i++)
						for(int j=0;j<m_buffer[i].size();j++)
							distances[i]+=(m_buffer[i][j]!=values[j]);

				}

				if(m_useSparseStorage == true)
				{
					// use sparse matrix impl
					vector<vector<float> > temp(1);
					temp[0] = distances;
					temp = m_populations[0]->MPI()->MPILayerReduceVariable(m_results);

					for(int i=0;i<temp.size();i++)
					{
						if(temp[0][i]<m_sparseThreshold)
						{
							// store (in final matrix)
						}
					}
				}
				else
				{
 					m_results.push_back(distances); // will be dense, triangular if checking against previous timesteps, rectangular if checking against specific dataset/training set
				}
			}
		}
	}
}

void AnalysisDistance::FinalizeSegmResults(vector<vector<float> > targets, int classStep, int dynStep)
{
	this->Finalize();

	if(m_results.size()>0)
	{
		vector<float> midRes = this->GetSegmOutputResult(targets,classStep,dynStep);
		if(m_results2.size() == 0)
			m_results2.push_back(midRes);
		else
		{
			for(int i=0;i<midRes.size();i++)
				m_results2[0].push_back(midRes[i]);
			//m_results2[m_results2.size()-1].insert(m_results2[m_results2.size()-1].end(),midRes);
		}

		m_results.clear(); // remove buffer
	}
}

// compares finalized results with output vector (of mixtures)
// 1st part, fully correct, yes or no
// 2nd part, how many correct was in there
// version 2: if there are mapped outputs, compare to these instead of index of closest item
// version 3: called from another function in order to do store analysis iteratively during simulation
vector<float> AnalysisDistance::GetSegmOutputResult(vector<vector<float> > targets, int classStep, int dynStep)
{
	int step = 0;
	int currentClass = 0;
	vector<float> result2;
	vector<float> result;

	cout<<"A";//cout<<"Targets.size() == "<<targets.size()<<", m_results.size() == "<<m_results.size()<<", classStep == "<<classStep<<", dynStep == "<<dynStep<<"...";
	cout.flush();

//	Storage storage;
	// debug test
	vector<vector<float> > vf(1);
	vector<float> f(m_specificCorrespondingClassIndexes.size());
	for(int i=0;i<m_specificCorrespondingClassIndexes.size();i++)
	{
		f[i]= (float)m_specificCorrespondingClassIndexes[i];
	}
	vf[0] = f;
	vector<float> indexesStorage;


	// end debug test

	for(int i=0;i<targets.size();i++)
	{
		vector<int> indexes;

		for(int j=dynStep-2;j<classStep;j+=dynStep)
		{
			vector<float> res = m_results[i*classStep+j];
			float lowVal = 1e8;
			int index = -1;
			for(int i=0;i<res.size();i++)
			{
				if(res[i]<=lowVal)
				{
					lowVal = res[i];
					index = i-1;// adding -1 to get right index
				}
			}

			if(index == -1) // should not happen
			{
				cout<<"Error:!-1";cout.flush();
			}

			// if mapped outputs exist, change to this
			if(this->m_specificCorrespondingClassIndexes.size()>0)
				index = m_specificCorrespondingClassIndexes[index];

			indexesStorage.push_back(index);

			for(int m=0;m<targets[i].size();m++)//currentClass].size();m++)
			{
			//	if(targets[i][m] == index)
//				{
					bool contains = false;
					for(int n=0;n<indexes.size();n++)
						if(indexes[n] == index)
							contains = true;
					if(contains == false)
						indexes.push_back(index);

				//	break;
				//}

			}

			step++;
		}

		int nrCorrect = 0;
		for(int m=0;m<targets[i].size();m++)
		{
			if(find(indexes.begin(),indexes.end(),(int)(targets[i][m])) != indexes.end())
			{
				nrCorrect++;
			}
		}
		
		if(nrCorrect == targets[i].size() && indexes.size() == targets[i].size())
			result.push_back(1.0);
		else if(nrCorrect == targets[i].size())
			result.push_back(0.0);//0.5);
		else
			result.push_back(0.0);

		result2.push_back((float)nrCorrect/(float)(targets[i].size()));

		//if(i>0 && i%classStep == 0)
			currentClass++;
	}

	//result = result / (float)targets.size();// only checkin if it is completely correct;

	for(int i=0;i<result2.size();i++)
		result.push_back(result2[i]);

//	cout<<"ok.\n";

	// debug test
	vf.push_back(indexesStorage);
//	storage.SaveDataFloatCSV("specificCorrespondingClassIndexes.csv",vf,Storage::Standard,"");

	return result;
}

void AnalysisDistance::Finalize()
{	
	if(m_useSparseStorage == true)
	{
		//// own file writing
		//// write sparse matrix to disk
	}
	else
	{
		/*if(m_results.size()==0)
		{
			cout<<"["<<m_populations[0]->network()->MPIGetNodeId()<<"]";cout.flush();
		}*/

		//cout<<this->m_populations[0]->network()->MPIGetNodeId()<<" ";cout.flush();
		m_results = m_populations[0]->MPI()->MPILayerReduceVariable(m_results); // current impl high memory usage
		//cout<<"{"<<this->m_populations[0]->network()->MPIGetNodeId()<<" ";cout.flush();
		/*if(m_measure == AnalysisDistance::Euclidean)
		{
			for(int i=0;i<m_results.size();i++)
			{

			for(int j=0;j<m_results[i].size();j++)
					m_results[i][j] = sqrt(m_results[i][j]);
			}
		}*/

		DoOwnFileWriting(); // if activated
	}
}