#include "NetworkDepression.h"

ConnectionModifierDepression::ConnectionModifierDepression()
{
	m_eventId = 16;
	m_strength = 0.1;
	m_hasStoredOrgWeights = false;
}


ConnectionModifierDepression::ConnectionModifierDepression(float strength)
{
	m_eventId = 16;
	m_strength = strength;
	m_hasStoredOrgWeights = false;
}

void ConnectionModifierDepression::SetConnection(Connection* c)
{
	m_connection = c;
}

void ConnectionModifierDepression::Modify()
{
	if(!IsOn()) return;

	float weight,diff;
	vector<long> postIds = m_connection->GetPostIds();
	vector<long> preIds;


	if(!m_hasStoredOrgWeights){
		//Store the original weights so we can reset, is done once.
		for(int i=0; i<postIds.size();i++){
			preIds = m_connection->GetPreIds(postIds[i]);
			m_orgWeights.push_back(vector<long> (preIds.size()));
			for(int j=0; j<preIds.size();j++){
				m_orgWeights[i][j] = m_connection->network()->GetWeight(preIds[j], postIds[i]);
			}
		}
		m_hasStoredOrgWeights = true;
	}

	for(int i=0; i<postIds.size();i++){
		preIds = m_connection->GetPreIds(postIds[i]);
		//Is the unit active?
		for(int j=0; j<preIds.size();j++){
			if(m_connection->network()->GetPreValue(preIds[j]) > 0){
				//depress
				weight = m_connection->network()->GetWeight(preIds[j], postIds[i]);
				diff = (weight-1);

				if(diff > 0){
					weight -= diff*m_strength;
				} else if(diff < 0){
					weight += -diff*m_strength;
				}

				m_connection->network()->SetWeight(weight,preIds[j], postIds[i]);

			} else {
				//reset weight
				m_connection->network()->SetWeight(m_orgWeights[i][j],preIds[j], postIds[i]);
			}
		}
	}
}


void ConnectionModifierDepression::Simulate(UnitModifier* e) {

}
