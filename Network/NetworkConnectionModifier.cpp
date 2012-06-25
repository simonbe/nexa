#include <math.h>
#include "NetworkConnectionModifier.h"
#include "NetworkUnits.h"
//#include "NetworkUnitModifier.h"



void ConnectionModifier::Initialize(Connection* connection)
{
	m_connection = connection;

	for(int i = 0;i<m_parentPopulationModifier.size();i++)
	{
		m_parentPopulationModifier[i]->AddChildConnectionModifier(this);
	}

	for(int i=0;i<m_parentConnectionModifier.size();i++)
	{
		m_parentConnectionModifier[i]->AddChildConnectionModifier(this);
	}
}

UnitModifier* ConnectionModifier::GetTransferFunction()
{
	return m_transferFunction;
}

ConnectionModifier::~ConnectionModifier()
{
	if(m_transferFunction!=NULL)
		delete m_transferFunction;
	
/*	for(int i=0;i<m_parentConnectionModifier.size();i++)
		delete m_parentConnectionModifier[i];

	for(int i=0;i<m_childConnectionModifier.size();i++)
		delete m_childConnectionModifier[i];

	for(int i=0;i<m_parentPopulationModifier.size();i++)
		delete m_parentPopulationModifier[i];
		*/
}